#ifndef _MEMIC_MODEL_H
#define _MEMIC_MODEL_H

#include "ResourceGradient.h"
#include "config/ArgManager.h"
#include "Evolve/World.h"


// Default values for plate dimensions are extracted from MEMIC plate stl 

EMP_BUILD_CONFIG( MemicConfig,
  GROUP(MAIN, "Global settings"),
  VALUE(SEED, int, -1, "Random number generator seed"),
  VALUE(TIME_STEPS, int, 1000, "Number of time steps to run for"),
  VALUE(PLATE_LENGTH, double, 10.0, "Length of plate in mm"),
  VALUE(PLATE_WIDTH, double, 6.0, "Width of plate in mm"),
  VALUE(PLATE_DEPTH, double, 1.45, "Depth of plate in mm"), 
  VALUE(CELL_DIAMETER, double, 20.0, "Cell length and width in microns"),
  // VALUE(WORLD_X, int, 100, "Width of world (in number of cells)"),
  // VALUE(WORLD_Y, int, 100, "Height of world (in number of cells)"),
  VALUE(INIT_POP_SIZE, int, 100, "Number of cells to seed population with"),

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
  VALUE(DIFFUSION_STEPS_PER_TIME_STEP, int, 10, "Rate at which diffusion is calculated relative to rest of model"),
  VALUE(OXYGEN_THRESHOLD, double, .1, "How much oxygen do cells need to survive?"),
  VALUE(KM, double, 0.01, "Michaelis-Menten kinetic parameter"),
);

struct Cell {
    double stemness = 1;
    int age = 0;
    int clade = 0;
    double hif1alpha = 0;

    Cell(double in_stemness = 0) : 
      stemness(in_stemness) {;}

    bool operator== (const Cell & other) const {
      return age == other.age && clade == other.clade;
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

  int WORLD_X;
  int WORLD_Y;
  int WORLD_Z;

  int next_clade = 1;

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

    WORLD_X = floor(PLATE_WIDTH / (CELL_DIAMETER/1000));
    WORLD_Y = floor(PLATE_LENGTH / (CELL_DIAMETER/1000));
    WORLD_Z = floor(PLATE_DEPTH / (CELL_DIAMETER/1000));

    if (oxygen) {
      oxygen->SetDiffusionCoefficient(OXYGEN_DIFFUSION_COEFFICIENT);
    }
  }

  int GetWorldX() {
    return WORLD_X;
  }

  int GetWorldY() {
    return WORLD_Y;
  }

  int GetWorldZ() {
    return WORLD_Z;
  }


  void InitPop() {
    // for (int cell_id = 0; cell_id < WORLD_X * WORLD_Y; cell_id++) {
    //   InjectAt(Cell(CELL_STATE::HEALTHY), cell_id);
    // }
    pop.resize(WORLD_X*WORLD_Y);
    int initial_spot = random_ptr->GetUInt(WORLD_Y*WORLD_X);
    InjectAt(Cell(), initial_spot);

    for (int cell_id = 0; cell_id < INIT_POP_SIZE; cell_id++) {
      int spot = random_ptr->GetUInt(WORLD_Y*WORLD_X);
      while (spot == initial_spot) {
        spot = random_ptr->GetUInt(WORLD_Y*WORLD_X);        
      }
      AddOrgAt(emp::NewPtr<Cell>(), spot, initial_spot);
    }

    RemoveOrgAt(initial_spot);
  }

  void InitOxygen() {
    for (int x = 0; x < WORLD_X; x++) {
      for (int y = 0; y < WORLD_Y; y++) {
        for (int z = 0; z < WORLD_Z; z++) {
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
      for (int x = 0; x < WORLD_X; x++) {
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

    SetPopStruct_Grid(WORLD_X, WORLD_Y, true);
    InitOxygen();
    InitPop();
  }

  void BasalOxygenConsumption() {
    for (int cell_id = 0; cell_id < (int)pop.size(); cell_id++) {
      if (IsOccupied(cell_id)) {
        int x = cell_id % WORLD_X;
        int y = cell_id / WORLD_X;
        double oxygen_loss_multiplier = oxygen->GetVal(x, y, 0);
        oxygen_loss_multiplier /= oxygen_loss_multiplier + KM;
        oxygen->DecNextVal(x, y, 0, BASAL_OXYGEN_CONSUMPTION * oxygen_loss_multiplier);
      }
    }
  }


  /// Determine if cell can divide (i.e. is space available). If yes, return
  /// id of cell that it can divide into. If not, return -1.
  int CanDivide(int cell_id) {

    emp::vector<int> open_spots;
    int x_coord = (cell_id % WORLD_X);
    int y_coord = (cell_id / WORLD_X);
    
    // Iterate over 9-cell neighborhood. Currently checks focal cell uneccesarily,
    // but that shouldn't cause problems because it will never show up as invasible.
    for (int x = std::max(0, x_coord-1); x < std::min(WORLD_X, x_coord + 2); x++) {
      for (int y = std::max(0, y_coord-1); y < std::min(WORLD_Y, y_coord + 2); y++) {
        int this_cell = y*WORLD_X + x;
        // Cells can be divided into if they are empty or if they are healthy and the
        // dividing cell is cancerous
        if (!IsOccupied(this_cell)) {
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


  void Quiesce(int cell_id) {
    // Quiescence - stick the cell back into the population in
    // the same spot but don't change anything else
    // std::cout << "Quieseing" << std::endl;
    pop[cell_id]->age++;
    emp::Ptr<Cell> cell = emp::NewPtr<Cell>(*pop[cell_id]);
    AddOrgAt(cell, emp::WorldPosition(cell_id,1), cell_id);
  }

  void RunStep() {

    std::cout << update << std::endl;

    for (int cell_id = 0; cell_id < WORLD_X * WORLD_Y; cell_id++) {

      if (!IsOccupied(cell_id)) {
        // Don't need to do anything for dead/empty cells
        continue;
      }

      int x = cell_id % WORLD_X;
      int y = cell_id / WORLD_X;

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
        // Cell divides
        oxygen->DecNextVal(x, y, 0, OXYGEN_CONSUMPTION_DIVISION);

        // Handle daughter cell in previously empty spot
        before_repro_sig.Trigger(cell_id);
        emp::Ptr<Cell> offspring = emp::NewPtr<Cell>(*pop[cell_id]);
        Mutate(offspring);
        offspring_ready_sig.Trigger(*offspring, cell_id);
        AddOrgAt(offspring, emp::WorldPosition(potential_offspring_cell, 1), cell_id);

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
  }

};




#endif