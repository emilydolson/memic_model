#ifndef _MEMIC_MODEL_H
#define _MEMIC_MODEL_H

#include "ResourceGradient.h"
#include "config/ArgManager.h"
#include "Evolve/World.h"


EMP_BUILD_CONFIG( MemicConfig,
  GROUP(MAIN, "Global settings"),
  VALUE(SEED, int, -1, "Random number generator seed"),
  VALUE(TIME_STEPS, int, 1000, "Number of time steps to run for"),
  VALUE(ASYMMETRIC_DIVISION_PROB, double, 0, "Probability of a change in stemness"),
  VALUE(MITOSIS_PROB, double, 0, "Probability of mitosis"),
  VALUE(OXYGEN_THRESHOLD, double, 0, "How much oxygen do cells need to survive?"),
  VALUE(HYPOXIA_DEATH_PROB, double, 0, "Probability of dieing, given hypoxic conditions"),
  VALUE(AGE_LIMIT, int, 0, "Age over which non-stem cells die"),
  VALUE(WORLD_X, int, 100, "Width of world (in number of cells)"),
  VALUE(WORLD_Y, int, 100, "Height of world (in number of cells)"),
);

enum class CELL_STATE {TUMOR=0, DEAD=2, HEALTHY=3};

struct Cell {
    double stemness;
    CELL_STATE state;
    int age = 0;

    Cell(CELL_STATE in_state, double in_stemness = 0) : 
      stemness(in_stemness), state(in_state) {;}
};

class HCAWorld : public emp::World<Cell> {
  private:
  int TIME_STEPS;
  double ASYMMETRIC_DIVISION_PROB;
  double MITOSIS_PROB;
  double OXYGEN_THRESHOLD;
  double HYPOXIA_DEATH_PROB;
  int AGE_LIMIT;
  int WORLD_X;
  int WORLD_Y;

  public:

  HCAWorld(emp::Random & r) : emp::World<Cell>(r) {;}

  void InitConfigs(MemicConfig & config) {
    TIME_STEPS = config.TIME_STEPS();
    ASYMMETRIC_DIVISION_PROB = config.ASYMMETRIC_DIVISION_PROB();
    MITOSIS_PROB = config.MITOSIS_PROB();
    OXYGEN_THRESHOLD = config.OXYGEN_THRESHOLD();
    HYPOXIA_DEATH_PROB = config.HYPOXIA_DEATH_PROB();
    AGE_LIMIT = config.AGE_LIMIT();
    WORLD_X = config.WORLD_X();
    WORLD_Y = config.WORLD_Y();
  }

  void InitPop() {
    for (int cell_id = 0; cell_id < WORLD_X * WORLD_Y; cell_id++) {
      InjectAt(Cell(CELL_STATE::HEALTHY), cell_id);
    }
  }

  void Setup(MemicConfig & config) {
    InitConfigs(config);

    // Setup systematics here

    SetPopStruct_Grid(WORLD_X, WORLD_Y, true);
    InitPop();
  }

  int CanDivide(int cell_id) {

    bool healthy = pop[cell_id]->state == CELL_STATE::HEALTHY;

    emp::vector<int> open_spots;
    int x_coord = (cell_id % WORLD_X);
    int y_coord = (cell_id / WORLD_X);
    for (int x = std::max(0, x_coord-1); x < std::min(WORLD_X, x_coord + 2); x++) {
      for (int y = std::max(0, y_coord-1); y < std::min(WORLD_Y, y_coord + 2); y++) {
        int this_cell = y*WORLD_X + x;
        if (!IsOccupied(this_cell) || pop[this_cell]->state == CELL_STATE::DEAD ||
            (pop[this_cell]->state == CELL_STATE::HEALTHY && !healthy)) {
          open_spots.push_back(this_cell);
        }
      }
    }
    
    if (open_spots.size() == 0) {
      return -1;
    }

    // If there are one or more available spaces, return a random spot
    return open_spots[random_ptr->GetUInt(0, open_spots.size())];
  }

  void RunStep() {
    for (int cell_id = 0; cell_id < WORLD_X * WORLD_Y; cell_id++) {

      // Query oxygen to test for hypoxia
      // If hypoxic, possibly die
      // Else, consume oxygen

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