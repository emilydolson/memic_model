#ifndef _MEMIC_MODEL_H
#define _MEMIC_MODEL_H

#include "ResourceGradient.h"
#include "config/ArgManager.h"
#include "tools/File.h"
#include "Evolve/World.h"
#include "tools/spatial_stats.h"

// Default values for plate dimensions are extracted from MEMIC plate stl 

EMP_BUILD_CONFIG( MemicConfig,
  GROUP(MAIN, "Global settings"),
  VALUE(SEED, int, -1, "Random number generator seed"),
  VALUE(TIME_STEPS, int, 1000, "Number of time steps to run for"),
  VALUE(PLATE_LENGTH, double, 10.0, "Length of plate in mm"),
  VALUE(PLATE_WIDTH, double, 6.0, "Width of plate in mm"),
  VALUE(PLATE_DEPTH, double, 1.45, "Depth of plate in mm"), 
  VALUE(CELL_DIAMETER, double, 20.0, "Cell length and width in microns"),
  VALUE(INIT_POP_SIZE, int, 100, "Number of cells to seed population with"),
  VALUE(DATA_RESOLUTION, int, 10, "How many updates between printing data?"),

  GROUP(CELL, "Cell settings"),
  VALUE(NEUTRAL_MUTATION_RATE, double, .05, "Probability of a neutral mutation (only relevant for phylogenetic signature)"),
  VALUE(ASYMMETRIC_DIVISION_PROB, double, 0, "Probability of a change in stemness"),
  VALUE(MITOSIS_PROB, double, .5, "Probability of mitosis"),
  VALUE(HYPOXIA_DEATH_PROB, double, .25, "Probability of dieing, given hypoxic conditions"),
  VALUE(AGE_LIMIT, int, 100, "Age over which non-stem cells die"),
  VALUE(BASAL_OXYGEN_CONSUMPTION, double, .00075, "Base oxygen consumption rate"),
  VALUE(OXYGEN_CONSUMPTION_DIVISION, double, .00075*5, "Amount of oxygen a cell consumes on division"),
  
  GROUP(OXYGEN, "Oxygen settings"),
  VALUE(INITIAL_OXYGEN_LEVEL, double, .5, "Initial oxygen level (will be placed in all cells)"),
  VALUE(OXYGEN_DIFFUSION_COEFFICIENT, double, .1, "Oxygen diffusion coefficient"),
  VALUE(DIFFUSION_STEPS_PER_TIME_STEP, int, 100, "Rate at which diffusion is calculated relative to rest of model"),
  VALUE(OXYGEN_THRESHOLD, double, .1, "How much oxygen do cells need to survive?"),
  VALUE(KM, double, 0.01, "Michaelis-Menten kinetic parameter"),

  GROUP(TREATMENT, "Treatment settings"),
  VALUE(DOSES, int, 0, "Number of radiation doses to apply"),
  VALUE(DOSE_SIZE, double, 2, "Dose size (Gy)"),
  VALUE(DOSE_TIME, int, -1, "When to apply radiation"),
  VALUE(K_OER, double, 3.28, "Effective OER constant"),  
  VALUE(OER_MIN, double, 1, "OER min constant"),  
  VALUE(OER_ALPHA_MAX, double, 1.75, "OER alpha max constant"),  
  VALUE(OER_BETA_MAX, double, 3.25, "OER alpha (? this is what the paper says but I feel like it's supposed to be beta) max constant"),  
);

struct Cell {
    double stemness = 1;
    int age = 0;
    int clade = 0;
    double hif1alpha = 0;
    bool marked_for_death = false;

    Cell(double in_stemness = 0) : 
      stemness(in_stemness) {;}

    bool operator== (const Cell & other) const {
      return age == other.age && clade == other.clade;
    }

    bool operator< (const Cell & other) const {
      return clade < other.clade;
    }

};

class HCAWorld : public emp::World<Cell> {
  protected:
  int TIME_STEPS;
  double NEUTRAL_MUTATION_RATE;
  double ASYMMETRIC_DIVISION_PROB;
  double OXYGEN_DIFFUSION_COEFFICIENT;
  double MITOSIS_PROB;
  double OXYGEN_THRESHOLD;
  double OXYGEN_CONSUMPTION_DIVISION;
  double HYPOXIA_DEATH_PROB;
  int AGE_LIMIT;
  int INIT_POP_SIZE;
  int DIFFUSION_STEPS_PER_TIME_STEP;
  double BASAL_OXYGEN_CONSUMPTION;
  double INITIAL_OXYGEN_LEVEL;
  double KM;

  double PLATE_LENGTH;
  double PLATE_WIDTH;
  double PLATE_DEPTH;
  double CELL_DIAMETER;

  int DOSES;
  double DOSE_SIZE;
  int DOSE_TIME;
  double K_OER;
  double OER_MIN;
  double OER_ALPHA_MAX;
  double OER_BETA_MAX;

  size_t WORLD_X;
  size_t WORLD_Y;
  size_t WORLD_Z;

  int next_clade = 1;

  emp::vector<emp::vector<double>> densities;
  emp::vector<emp::vector<double>> diversities;

  std::function<double()> colless_fun = [this](){return GetSystematics()->CollessLikeIndex();};
  std::function<double()> sackin_fun = [this](){return GetSystematics()->SackinIndex();};

  public:
  emp::Ptr<ResourceGradient> oxygen;

  HCAWorld(emp::Random & r) : emp::World<Cell>(r), oxygen(nullptr) {;}
  HCAWorld() {;}

  ~HCAWorld() {
    if (oxygen) {
      oxygen.Delete();
    }
  }

  void InitConfigs(MemicConfig & config) {
    TIME_STEPS = config.TIME_STEPS();
    NEUTRAL_MUTATION_RATE = config.NEUTRAL_MUTATION_RATE();
    ASYMMETRIC_DIVISION_PROB = config.ASYMMETRIC_DIVISION_PROB();
    MITOSIS_PROB = config.MITOSIS_PROB();
    OXYGEN_DIFFUSION_COEFFICIENT = config.OXYGEN_DIFFUSION_COEFFICIENT();
    OXYGEN_THRESHOLD = config.OXYGEN_THRESHOLD();
    OXYGEN_CONSUMPTION_DIVISION = config.OXYGEN_CONSUMPTION_DIVISION();
    HYPOXIA_DEATH_PROB = config.HYPOXIA_DEATH_PROB();
    AGE_LIMIT = config.AGE_LIMIT();
    INITIAL_OXYGEN_LEVEL = config.INITIAL_OXYGEN_LEVEL();
    DIFFUSION_STEPS_PER_TIME_STEP = config.DIFFUSION_STEPS_PER_TIME_STEP();
    BASAL_OXYGEN_CONSUMPTION = config.BASAL_OXYGEN_CONSUMPTION();
    KM = config.KM();
    INIT_POP_SIZE = config.INIT_POP_SIZE();
    PLATE_LENGTH = config.PLATE_LENGTH();
    PLATE_WIDTH = config.PLATE_WIDTH();
    PLATE_DEPTH = config.PLATE_DEPTH();
    CELL_DIAMETER = config.CELL_DIAMETER();

    DOSES = config.DOSES();
    DOSE_SIZE = config.DOSE_SIZE();
    DOSE_TIME = config.DOSE_TIME();
    K_OER = config.K_OER();
    OER_MIN = config.OER_MIN();
    OER_ALPHA_MAX = config.OER_ALPHA_MAX();
    OER_BETA_MAX = config.OER_BETA_MAX();

    WORLD_X = (size_t)floor(PLATE_WIDTH / (CELL_DIAMETER/1000));
    WORLD_Y = (size_t)floor(PLATE_LENGTH / (CELL_DIAMETER/1000));
    WORLD_Z = (size_t)floor(PLATE_DEPTH / (CELL_DIAMETER/1000));

    if (oxygen) {
      oxygen->SetDiffusionCoefficient(OXYGEN_DIFFUSION_COEFFICIENT);
    }
  }

  size_t GetWorldX() {
    return WORLD_X;
  }

  size_t GetWorldY() {
    return WORLD_Y;
  }

  size_t GetWorldZ() {
    return WORLD_Z;
  }


  void InitPop() {
    // for (int cell_id = 0; cell_id < WORLD_X * WORLD_Y; cell_id++) {
    //   InjectAt(Cell(CELL_STATE::HEALTHY), cell_id);
    // }
    pop.resize(WORLD_X*WORLD_Y);
    size_t initial_spot = random_ptr->GetUInt(WORLD_Y*WORLD_X);
    InjectAt(Cell(), initial_spot);

    for (size_t cell_id = 0; cell_id < (size_t)INIT_POP_SIZE; cell_id++) {
      size_t spot = random_ptr->GetUInt(WORLD_Y*WORLD_X);
      while (spot == initial_spot) {
        spot = random_ptr->GetUInt(WORLD_Y*WORLD_X);        
      }
      AddOrgAt(emp::NewPtr<Cell>(), spot, initial_spot);
    }

    RemoveOrgAt(initial_spot);
  }

  void InitOxygen() {
    for (size_t x = 0; x < WORLD_X; x++) {
      for (size_t y = 0; y < WORLD_Y; y++) {
        for (size_t z = 0; z < WORLD_Z; z++) {
          oxygen->SetVal(x, y, z, INITIAL_OXYGEN_LEVEL);
        }
      }
    }
  }

  ResourceGradient& GetOxygen() {
    return *oxygen;
  }

  void UpdateOxygen() {
      BasalOxygenConsumption();
      oxygen->Diffuse();
      oxygen->Update();

      // Oxygen inflow along edge
      for (size_t x = 0; x < WORLD_X; x++) {
        oxygen->SetVal(x, 0, WORLD_Z-1, 1);
      }
  }

  void Reset(MemicConfig & config, bool web = false) {
    Clear();
    if (oxygen) {
      oxygen.Delete();
      oxygen = nullptr;
    }
    Setup(config, web);    
  }

  void Setup(MemicConfig & config, bool web = false) {
    InitConfigs(config);
    oxygen.New(WORLD_X, WORLD_Y, WORLD_Z);
    oxygen->SetDiffusionCoefficient(OXYGEN_DIFFUSION_COEFFICIENT);

    if (!web) { // Web version needs to do diffusion separately to visualize
      OnUpdate([this](int ud){
        for (int i = 0; i < DIFFUSION_STEPS_PER_TIME_STEP; i++) {
          UpdateOxygen();
        }
      });
    }

    emp::Ptr<emp::Systematics<Cell, int> > sys;
    sys.New([](const Cell & c){return c.clade;});
    AddSystematics(sys);

    // SetupFitnessFile().SetTimingRepeat(config.DATA_RESOLUTION());
    SetupSystematicsFile().SetTimingRepeat(config.DATA_RESOLUTION());
    SetupPopulationFile().SetTimingRepeat(config.DATA_RESOLUTION());

    emp::DataFile & phylodiversity_file = SetupFile("phylodiversity.csv");
    sys->AddEvolutionaryDistinctivenessDataNode();
    sys->AddPairwiseDistanceDataNode();
    sys->AddPhylogeneticDiversityDataNode();

    phylodiversity_file.AddVar(update, "generation", "Generation");
    phylodiversity_file.AddStats(*sys->GetDataNode("evolutionary_distinctiveness") , "evolutionary_distinctiveness", "evolutionary distinctiveness for a single update", true, true);
    phylodiversity_file.AddStats(*sys->GetDataNode("pairwise_distance"), "pairwise_distance", "pairwise distance for a single update", true, true);
    phylodiversity_file.AddCurrent(*sys->GetDataNode("phylogenetic_diversity"), "current_phylogenetic_diversity", "current phylogenetic diversity", true, true);

    phylodiversity_file.AddFun(colless_fun, "colless_index", "current colless index");  
    phylodiversity_file.AddFun(sackin_fun, "sackin_index", "current sackin index");  

    phylodiversity_file.PrintHeaderKeys();
    phylodiversity_file.SetTimingRepeat(config.DATA_RESOLUTION());

    // emp::DataFile & phylodiversity_file = SetupFile("phylodiversity.csv");

    SetPopStruct_Grid(WORLD_X, WORLD_Y, true);
    InitOxygen();
    InitPop();

    SetSynchronousSystematics(true);
  }

  void BasalOxygenConsumption() {
    for (size_t cell_id = 0; cell_id < pop.size(); cell_id++) {
      if (IsOccupied(cell_id)) {
        size_t x = cell_id % WORLD_X;
        size_t y = cell_id / WORLD_X;
        double oxygen_loss_multiplier = oxygen->GetVal(x, y, 0);
        oxygen_loss_multiplier /= oxygen_loss_multiplier + KM;
        oxygen->DecNextVal(x, y, 0, BASAL_OXYGEN_CONSUMPTION * oxygen_loss_multiplier);
      }
    }
  }


  /// Determine if cell can divide (i.e. is space available). If yes, return
  /// id of cell that it can divide into. If not, return -1.
  int CanDivide(size_t cell_id) {
    emp::vector<int> open_spots;
    int x_coord = (int)(cell_id % WORLD_X);
    int y_coord = (int)(cell_id / WORLD_X);
    
    // Iterate over 9-cell neighborhood. Currently checks focal cell uneccesarily,
    // but that shouldn't cause problems because it will never show up as invasible.
    for (int x = std::max(0, x_coord-1); x < std::min((int)WORLD_X, x_coord + 2); x++) {
      for (int y = std::max(0, y_coord-1); y < std::min((int)WORLD_Y, y_coord + 2); y++) {
        int this_cell = y*(int)WORLD_X + x;
        // Cells can be divided into if they are empty or if they are healthy and the
        // dividing cell is cancerous
        if (!IsOccupied((size_t)this_cell)) {
          open_spots.push_back(this_cell);
        }
      }
    }
    
    // -1 is a sentinel value indicating no spots are available
    if (open_spots.size() == 0) {
      return -1;
    }

    // If there are one or more available spaces, return a random spot
    return open_spots[random_ptr->GetUInt(0, open_spots.size())];
  }

  int Mutate(emp::Ptr<Cell> c){
    c->age = 0;
    
    if (random_ptr->P(NEUTRAL_MUTATION_RATE)) {
      c->clade = next_clade;
      next_clade++;
      return 1;
    }

    return 0;
  }


  void Quiesce(size_t cell_id) {
    // Quiescence - stick the cell back into the population in
    // the same spot but don't change anything else
    // std::cout << "Quieseing" << std::endl;
    pop[cell_id]->age++;
    if (pop[cell_id]->age < AGE_LIMIT) {
      emp::Ptr<Cell> cell = emp::NewPtr<Cell>(*pop[cell_id]);
      AddOrgAt(cell, emp::WorldPosition(cell_id,1), cell_id);      
    }
  }

  void RunStep() {
    std::cout << update << std::endl;

    if ((int)update == DOSE_TIME) {
      ApplyRadiation(DOSES, DOSE_SIZE);
    }

    for (size_t cell_id = 0; cell_id < WORLD_X * WORLD_Y; cell_id++) {
      if (!IsOccupied(cell_id)) {
        // Don't need to do anything for dead/empty cells
        continue;
      }

      size_t x = cell_id % WORLD_X;
      size_t y = cell_id / WORLD_X;

      // Query oxygen to test for hypoxia
      if (oxygen->GetVal(x, y, 0) < OXYGEN_THRESHOLD) {
        // If hypoxic, the hif1-alpha surpressor gets turned off
        // causing hif1-alpha to accumulate
        pop[cell_id]->hif1alpha = 1;

        // If hypoxic, die with specified probability
        if (!random_ptr->P(HYPOXIA_DEATH_PROB)) {
          Quiesce(cell_id); // Cell survives to next generation
        }
        continue; // Division not allowed under hypoxia 
        // TODO: Consider replacing this with a function relating
        // division probability to oxygen availability, as in 
        // Grimes et al 2018
      } else {
        // If not hypoxic, the hif1-alpha surpressor is on
        // causing hif1-alpha to not accumulate
        pop[cell_id]->hif1alpha = 0;
      }

      // Check for space for division
      int potential_offspring_cell = CanDivide(cell_id);

      // If space, divide
      if (potential_offspring_cell != -1 && random_ptr->P(MITOSIS_PROB)) {
        // Check if cell needs to die
        if (pop[cell_id]->marked_for_death) {
          continue;
        }
        
        // Cell divides
        oxygen->DecNextVal(x, y, 0, OXYGEN_CONSUMPTION_DIVISION);

        // Handle daughter cell in previously empty spot
        before_repro_sig.Trigger(cell_id);
        emp::Ptr<Cell> offspring = emp::NewPtr<Cell>(*pop[cell_id]);
        Mutate(offspring);
        offspring_ready_sig.Trigger(*offspring, cell_id);
        AddOrgAt(offspring, emp::WorldPosition((size_t)potential_offspring_cell, 1), cell_id);

        // Handle daughter cell in current location
        before_repro_sig.Trigger(cell_id);
        offspring = emp::NewPtr<Cell>(*pop[cell_id]);
        Mutate(offspring);
        offspring_ready_sig.Trigger(*offspring, cell_id);
        AddOrgAt(offspring, emp::WorldPosition(cell_id,1), cell_id);
        // std::cout << "Mutated: " << offspring->clade << std::endl;
      } else {        
        Quiesce(cell_id);
      }
    }

    Update();
  }

  void Run() {
      for (int u = 0; u <= TIME_STEPS; u++) {
          RunStep();
      }
      PrintOxygenGrid("oxygen.csv");
      systematics[0].DynamicCast<emp::Systematics<Cell, int>>()->Snapshot("memic_phylo.csv");  
      densities = emp::GridDensity(*this);
      std::cout << emp::to_string(densities) << std::endl;

  }

  /* n is the number of doses of radiation, d dose size in Gy*/
  double SurvivingFraction(double n, double d, double c) {
    double alpha = OER_ALPHA_MAX/((((OER_ALPHA_MAX - OER_MIN)*K_OER)/(c + K_OER)) + OER_MIN);
    double beta = OER_BETA_MAX/emp::Pow(((((OER_BETA_MAX - OER_MIN)*K_OER)/(c + K_OER)) + OER_MIN), 2);
    return exp(-n*(alpha*d + beta*emp::Pow(d,2)));
  }

  /* n is the number of doses of radiation, d dose size in Gy*/
  void ApplyRadiation(double n, double d) {
    for (size_t cell_id = 0; cell_id < WORLD_X * WORLD_Y; cell_id++) {
      if (!IsOccupied(cell_id)) {
        // Don't need to do anything for dead/empty cells
        continue;
      }

      size_t x = cell_id % WORLD_X;
      size_t y = cell_id / WORLD_X;
      double c = oxygen->GetVal(x, y, 0);


      if (!random_ptr->P(SurvivingFraction(n,d,c))) {
        // TODO: Figure out best way to kill cells
        pop[cell_id]->marked_for_death = true;
      } 

    }
  }

  void PrintOxygenGrid(const std::string & filename) const {
    std::cout << "print o2" <<std::endl;
    std::ofstream oxygen_file(filename);

    for (size_t cell_id = 0; cell_id < WORLD_X * WORLD_Y; cell_id++) {
      size_t x = cell_id % WORLD_X;
      size_t y = cell_id / WORLD_X;
      if (x % WORLD_X > 0 ) {
        oxygen_file << ", "; // Don't add comma at beginning of line
      }

      oxygen_file << emp::to_string(oxygen->GetVal(x, y, 0));

      if (x % WORLD_X == WORLD_X - 1 ) {
        oxygen_file << "\n"; // We're at the end of a row
      }
    }

    oxygen_file.close();
  }
};

#endif