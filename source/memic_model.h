#ifndef _MEMIC_MODEL_H
#define _MEMIC_MODEL_H

#include "ResourceGradient.h"
#include "config/ArgManager.h"
#include "Evolve/World.h"


EMP_BUILD_CONFIG( MemicConfig,
  GROUP(MAIN, "Global settings"),
  VALUE(SEED, int, -1, "Random number generator seed"),
  VALUE(TIME_STEPS, int, 1000, "Number of time steps to run for"),
  VALUE(WORLD_X, int, 100, "Width of world (in number of cells)"),
  VALUE(WORLD_Y, int, 100, "Height of world (in number of cells)"),

  GROUP(CELL, "Cell settings"),
  VALUE(NEUTRAL_MUTATION_RATE, double, .01, "Probability of a neutral mutation (only relevant for phylogenetic signature)"),
  VALUE(ASYMMETRIC_DIVISION_PROB, double, 0, "Probability of a change in stemness"),
  VALUE(MITOSIS_PROB, double, 0, "Probability of mitosis"),
  VALUE(HYPOXIA_DEATH_PROB, double, 0, "Probability of dieing, given hypoxic conditions"),
  VALUE(AGE_LIMIT, int, 100, "Age over which non-stem cells die"),
  VALUE(BASAL_OXYGEN_CONSUMPTION_HEALTHY, double, .000375, "Base oxygen consumption rate for healthy cells"),
  VALUE(BASAL_OXYGEN_CONSUMPTION_TUMOR, double, .00075, "Base oxygen consumption rate for tumor cells"),

  GROUP(OXYGEN, "Oxygen settings"),
  VALUE(OXYGEN_DIFFUSION_COEFFICIENT, double, 0, "Oxygen diffusion coefficient"),
  VALUE(DIFFUSION_STEPS_PER_TIME_STEP, int, 10, "Rate at which diffusion is calculated relative to rest of model"),
  VALUE(OXYGEN_THRESHOLD, double, 0, "How much oxygen do cells need to survive?"),
  VALUE(OXYGEN_CONSUMPTION, double, 0, "Amount of oxygen a cell consumes per time point"),
  VALUE(KM, double, 0.01, "Michaelis-Menten kinetic parameter"),
);

enum class CELL_STATE {TUMOR=0, HEALTHY=1};

struct Cell {
    double stemness;
    CELL_STATE state;
    int age = 0;
    int clade = -1;

    Cell(CELL_STATE in_state, double in_stemness = 0) : 
      stemness(in_stemness), state(in_state) {;}
};

class HCAWorld : public emp::World<Cell> {
  private:
  int TIME_STEPS;
  double NEUTRAL_MUTATION_RATE;
  double ASYMMETRIC_DIVISION_PROB;
  double OXYGEN_DIFFUSION_COEFFICIENT;
  double MITOSIS_PROB;
  double OXYGEN_THRESHOLD;
  double OXYGEN_CONSUMPTION;
  double HYPOXIA_DEATH_PROB;
  int AGE_LIMIT;
  int WORLD_X;
  int WORLD_Y;
  int DIFFUSION_STEPS_PER_TIME_STEP;
  double BASAL_OXYGEN_CONSUMPTION_HEALTHY;
  double BASAL_OXYGEN_CONSUMPTION_TUMOR;
  double KM;

  int next_clade = 0;

  public:

  emp::Ptr<ResourceGradient> oxygen;

  HCAWorld(emp::Random & r) : emp::World<Cell>(r) {;}

  void InitConfigs(MemicConfig & config) {
    TIME_STEPS = config.TIME_STEPS();
    NEUTRAL_MUTATION_RATE = config.NEUTRAL_MUTATION_RATE();
    ASYMMETRIC_DIVISION_PROB = config.ASYMMETRIC_DIVISION_PROB();
    MITOSIS_PROB = config.MITOSIS_PROB();
    OXYGEN_DIFFUSION_COEFFICIENT = config.OXYGEN_DIFFUSION_COEFFICIENT();
    OXYGEN_THRESHOLD = config.OXYGEN_THRESHOLD();
    OXYGEN_CONSUMPTION = config.OXYGEN_CONSUMPTION();
    HYPOXIA_DEATH_PROB = config.HYPOXIA_DEATH_PROB();
    AGE_LIMIT = config.AGE_LIMIT();
    WORLD_X = config.WORLD_X();
    WORLD_Y = config.WORLD_Y();
    DIFFUSION_STEPS_PER_TIME_STEP = config.DIFFUSION_STEPS_PER_TIME_STEP();
    BASAL_OXYGEN_CONSUMPTION_HEALTHY = config.BASAL_OXYGEN_CONSUMPTION_HEALTHY();
    BASAL_OXYGEN_CONSUMPTION_TUMOR = config.BASAL_OXYGEN_CONSUMPTION_TUMOR();
    KM = config.KM();
  }

  void InitPop() {
    for (int cell_id = 0; cell_id < WORLD_X * WORLD_Y; cell_id++) {
      InjectAt(Cell(CELL_STATE::HEALTHY), cell_id);
    }
  }

  void Setup(MemicConfig & config) {
    InitConfigs(config);
    oxygen.New(WORLD_X, WORLD_Y);
    oxygen->SetDiffusionCoefficient(OXYGEN_DIFFUSION_COEFFICIENT);
    OnUpdate([this](int ud){
      for (int i = 0; i < DIFFUSION_STEPS_PER_TIME_STEP; i++) {
        BasalOxygenConsumption();
        oxygen->Diffuse();
        oxygen->Update();
      }
    });

    // For now, mutations are just a way for there to be
    // discernable phylogenetic signal.
    SetMutFun([this](Cell & c, emp::Random & rnd){
      if (c.state == CELL_STATE::TUMOR) {
        if (rnd.P(NEUTRAL_MUTATION_RATE)) {
          c.clade = next_clade;
          next_clade++;
          return 1;
        }
      }
      return 0;
    });

    SetAutoMutate();

    SetPopStruct_Grid(WORLD_X, WORLD_Y, true);
    InitPop();
  }

  void BasalOxygenConsumption() {
    for (int cell_id = 0; cell_id < (int)pop.size(); cell_id++) {
      if (IsOccupied(cell_id)) {
        int x = cell_id % WORLD_X;
        int y = cell_id / WORLD_X;
        double oxygen_loss_multiplier = oxygen->GetVal(x, y);
        oxygen_loss_multiplier /= oxygen_loss_multiplier + KM;

        switch (pop[cell_id]->state) {
          case CELL_STATE::HEALTHY:
            oxygen->DecNextVal(x, y, BASAL_OXYGEN_CONSUMPTION_HEALTHY * oxygen_loss_multiplier);
            break;

          case CELL_STATE::TUMOR:
            oxygen->DecNextVal(x, y, BASAL_OXYGEN_CONSUMPTION_TUMOR * oxygen_loss_multiplier);
            break;

          default:
            emp_assert(false, "INVALID CELL STATE");
            break;
        }
      }
    }
  }


  /// Determine if cell can divide (i.e. is space available). If yes, return
  /// id of cell that it can divide into. If not, return -1.
  int CanDivide(int cell_id) {

    bool healthy = pop[cell_id]->state == CELL_STATE::HEALTHY;

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
        if (!IsOccupied(this_cell) ||
            (pop[this_cell]->state == CELL_STATE::HEALTHY && !healthy)) {
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
      if (oxygen->GetVal(x, y) < OXYGEN_THRESHOLD) {
        // If hypoxic, die with specified probability
        if (random_ptr->P(HYPOXIA_DEATH_PROB)) {
          continue; // Continuing without reproducing means death on next update
        }
      }

      // Else, consume oxygen
      oxygen->DecNextVal(x, y, OXYGEN_CONSUMPTION);

      // Check for space for division
      int potential_offspring_cell = CanDivide(cell_id);

      // If space, divide
      if (potential_offspring_cell != -1) {
        // Cell divides

        // Handle daughter cell in previously empty spot
        before_repro_sig.Trigger(cell_id);
        emp::Ptr<Cell> offspring = emp::NewPtr<Cell>(*pop[cell_id]);
        offspring_ready_sig.Trigger(*offspring, cell_id);
        AddOrgAt(offspring, potential_offspring_cell, cell_id);

        // Handle daughter cell in current location
        before_repro_sig.Trigger(cell_id);
        offspring = emp::NewPtr<Cell>(*pop[cell_id]);
        offspring_ready_sig.Trigger(*offspring, cell_id);
        AddOrgAt(offspring, cell_id, cell_id);
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