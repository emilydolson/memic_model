#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"
#include "../source/ResourceGradient.h"
#include "../source/memic_model.h"

int y_len = 50;
int x_len = 100;
int z_len = 1;
MemicConfig config;
HCAWorld world;

TEST_CASE("Test constructor", "[oxygen_gradient]") {
    ResourceGradient r(x_len, y_len, z_len);
    
    // Make sure r got constructed correctly
    for (int y = 0; y < y_len; y++) {
        for (int x = 0; x < x_len; x++) {
            for (int z = 0; z < z_len; z++) {
                CHECK(r.GetVal(x, y,z) == 0);
            }
        }
    }

    // Test constructing from grid
    emp::vector<emp::vector<emp::vector<double>>> grid;
    int grid_y = 5;
    int grid_x = 3;
    int grid_z = 1;
    grid.resize(grid_z);
    for (int z = 0; z < grid_z; z++) {
        grid[z].resize(grid_y);
        for (int y = 0; y < grid_y; y++) {
            grid[z][y].resize(grid_x);
            for (int x = 0; x < grid_x; x++) {
                grid[z][y][x] = x*x+y;
            }
        }
    }

    ResourceGradient r2(grid);
    for (int y = 0; y < grid_y; y++) {
        for (int x = 0; x < grid_x; x++) {
            CHECK(r2.GetVal(x, y, 0) == x*x+y);
        }
    }

}

TEST_CASE("Test standard diffusion", "[oxygen_gradient]") {
    ResourceGradient r(x_len, y_len);

    // If there is no oxygen, there should still be no oxygen
    // after diffusion happens
    r.SetDiffusionCoefficient(.01);
    r.Diffuse();
    for (int y = 0; y < y_len; y++) {
        for (int x = 0; x < x_len; x++) {
            CHECK(r.GetVal(x, y, 0) == 0);
            CHECK(r.GetNeighborOxygen(x,y, 0) == 0);
        }
    }
    
    // Add some oxygen (and make sure that it works)
    r.SetVal(0,0,0,1);
    CHECK(r.GetVal(0,0,0) == 1);

    // Test GetNeighborOxygen
    // All the oxygen is in 0,0, not its neighbors,
    // But Dirichlet boundaries cause us to assume the cell
    // has four neighbors off-screen with equivalent oxygen
    CHECK(r.GetNeighborOxygen(0,0,0) == 4);
    CHECK(r.GetNeighborOxygen(1,0,0) == 1);
    CHECK(r.GetNeighborOxygen(0,1,0) == 1);
    CHECK(r.GetNeighborOxygen(1,1,0) == 0);

    // Not toroidal
    CHECK(r.GetNeighborOxygen(99,49,0) == 0);
    CHECK(r.GetNeighborOxygen(99,0,0) == 0);
    CHECK(r.GetNeighborOxygen(0,49,0) == 0);

    // Now do actual diffusion
    r.Diffuse();
    CHECK(r.GetNextVal(0,0,0) == Approx(0.98));
    CHECK(r.GetNextVal(0,1,0) == Approx(0.01));
    CHECK(r.GetNextVal(1,0,0) == Approx(0.01));
    CHECK(r.GetNextVal(1,1,0) == 0);

    // Not toroidal
    CHECK(r.GetNextVal(99,49,0) == 0);
    CHECK(r.GetNextVal(99,0,0) == 0);
    CHECK(r.GetNextVal(0,49,0) == 0);

}

TEST_CASE("Test turning off diffusion", "[oxygen_gradient]") {
    ResourceGradient r(x_len, y_len);

    // If there is no oxygen, there should still be no oxygen
    // after diffusion happens
    r.SetDiffusionCoefficient(0);
    r.Diffuse();
    for (int y = 0; y < y_len; y++) {
        for (int x = 0; x < x_len; x++) {
            CHECK(r.GetVal(x, y, 0) == 0);
            CHECK(r.GetNeighborOxygen(x,y,0) == 0);
        }
    }
    
    // Add some oxygen (and make sure that it works)
    r.SetVal(0,0,0,1);
    CHECK(r.GetVal(0,0,0) == 1);

    // Test GetNeighborOxygen
    // All the oxygen is in 0,0, not its neighbors
    CHECK(r.GetNeighborOxygen(0,0,0) == 4);
    CHECK(r.GetNeighborOxygen(1,0,0) == 1);
    CHECK(r.GetNeighborOxygen(0,1,0) == 1);
    CHECK(r.GetNeighborOxygen(1,1,0) == 0);

    // Not toroidal
    CHECK(r.GetNeighborOxygen(99,49,0) == 0);
    CHECK(r.GetNeighborOxygen(99,0,0) == 0);
    CHECK(r.GetNeighborOxygen(0,49,0) == 0);

    // Now do actual diffusion
    r.Diffuse();
    CHECK(r.GetNextVal(0,0,0) == Approx(1));
    CHECK(r.GetNextVal(0,1,0) == Approx(0));
    CHECK(r.GetNextVal(1,0,0) == Approx(0));
    CHECK(r.GetNextVal(1,1,0) == 0);

    // Not toroidal
    CHECK(r.GetNextVal(99,49,0) == 0);
    CHECK(r.GetNextVal(99,0,0) == 0);
    CHECK(r.GetNextVal(0,49,0) == 0);

}

TEST_CASE("Test toroidal diffusion", "[oxygen_gradient]") {
    ResourceGradient r(x_len, y_len);

    // If there is no oxygen, there should still be no oxygen
    // after diffusion happens
    r.SetDiffusionCoefficient(.01);
    r.SetToroidal(true);
    r.Diffuse();
    for (int y = 0; y < y_len; y++) {
        for (int x = 0; x < x_len; x++) {
            CHECK(r.GetVal(x, y,0) == 0);
            CHECK(r.GetNeighborOxygen(x,y,0) == 0);
        }
    }
    
    // Add some oxygen (and make sure that it works)
    r.SetVal(0,0,0,1);
    CHECK(r.GetVal(0,0,0) == 1);

    // Test GetNeighborOxygen
    // All the oxygen is in 0,0, not its neighbors,
    // but because this is toroidal it's its own
    // neighbor above and below
    CHECK(r.GetNeighborOxygen(0,0,0) == 2.0);
    CHECK(r.GetNeighborOxygen(1,0,0) == 1);
    CHECK(r.GetNeighborOxygen(0,1,0) == 1);
    CHECK(r.GetNeighborOxygen(1,1,0) == 0);

    // Toroidal
    CHECK(r.GetNeighborOxygen(99,49,0) == 0);
    CHECK(r.GetNeighborOxygen(99,0,0) == 1);
    CHECK(r.GetNeighborOxygen(0,49,0) == 1);

    // Now do actual diffusion
    r.Diffuse();
    CHECK(r.GetNextVal(0,0,0) == Approx(0.96));
    CHECK(r.GetNextVal(0,1,0) == Approx(0.01));
    CHECK(r.GetNextVal(1,0,0) == Approx(0.01));
    CHECK(r.GetNextVal(1,1,0) == 0);

    // Toroidal
    CHECK(r.GetNextVal(99,49,0) == 0);
    CHECK(r.GetNextVal(99,0,0) == 0.01);
    CHECK(r.GetNextVal(0,49,0) == 0.01);

}

TEST_CASE("Test Z axis", "[oxygen_gradient]") {
    ResourceGradient r(x_len, y_len, 10);

    // Add some oxygen (and make sure that it works)
    r.SetVal(1,2,3,2);
    CHECK(r.GetVal(1,2,3) == 2);

    CHECK(r.GetNeighborOxygen(1,2,3) == 0);
    CHECK(r.GetNeighborOxygen(1,2,2) == 2);
    CHECK(r.GetNeighborOxygen(1,2,4) == 2);
    CHECK(r.GetNeighborOxygen(1,1,3) == 2);
    CHECK(r.GetNeighborOxygen(1,3,3) == 2);
    CHECK(r.GetNeighborOxygen(2,2,3) == 2);
    CHECK(r.GetNeighborOxygen(0,2,3) == 2);

    r.SetToroidal(true);

    CHECK(r.GetNeighborOxygen(1,2,3) == 0);
    CHECK(r.GetNeighborOxygen(1,2,2) == 2);
    CHECK(r.GetNeighborOxygen(1,2,4) == 2);
    CHECK(r.GetNeighborOxygen(1,1,3) == 2);
    CHECK(r.GetNeighborOxygen(1,3,3) == 2);
    CHECK(r.GetNeighborOxygen(2,2,3) == 2);
    CHECK(r.GetNeighborOxygen(0,2,3) == 2);

}

TEST_CASE("Test updating gradient", "[oxygen_gradient]") {
    ResourceGradient r(x_len, y_len, 10);

    r.SetVal(10,10,0,5);
    r.SetVal(1,1,0,5);
    CHECK(r.GetVal(10,10,0) == 5);
    CHECK(r.GetVal(1,1,0) == 5);

    r.DecVal(10,10,0, 1);
    r.SetNextVal(8,6,0,3);
    r.SetNextVal(10,10,0,20);
    r.DecNextVal(10,10,0,2);
    r.SetNextVal(5,5,5,-4);
    CHECK(r.GetVal(10,10,0) == 4);
    CHECK(r.GetNextVal(8,6,0) == 3);
    CHECK(r.GetNextVal(10,10,0) == 18);
    CHECK(r.GetNextVal(5,5,5) == -4);

    r.Update();
    CHECK(r.GetVal(8,6,0) == 3);
    CHECK(r.GetVal(10,10,0) == 18);
    CHECK(r.GetVal(1,1,0) == 0);
    CHECK(r.GetVal(5,5,5) == 0); // Negative numbers should be zeroed out
}

TEST_CASE("Test HCAWorld", "[full_model]") {

    // Test destructor
    emp::Ptr<HCAWorld> world_ptr;
    emp::Random r;
    world_ptr.New(r);
    world_ptr.Delete();

    config.TIME_STEPS(5);
    config.CELL_DIAMETER(200);
    config.NEUTRAL_MUTATION_RATE(.75);
    world.Setup(config);

    CHECK(world.GetWorldX() == 30);
    CHECK(world.GetWorldY() == 50);
    CHECK(world.GetWorldZ() == 7);

    world.BasalOxygenConsumption();

    size_t example_id = 0;
    for (size_t cell_id = 0; cell_id < world.GetSize(); cell_id++) {
        if (world.IsOccupied(cell_id)) {
            example_id = cell_id;
            int age = world.GetOrg(cell_id).age;
            CHECK(!world.IsOccupied(emp::WorldPosition(cell_id, 1)));
            world.Quiesce(cell_id);
            CHECK(age + 1 == world.GetOrg(cell_id).age);
            CHECK(world.GetOrg(cell_id) == world.GetNextOrg(cell_id));

            Cell& cell_org = world.GetOrg(cell_id);
            int clade_before = cell_org.clade;            

            bool mutated = false;
            while(!mutated) {
                mutated = world.Mutate(&cell_org);
            }

            CHECK(clade_before != cell_org.clade);
            CHECK(cell_org.age == 0);

            CHECK(world.GetOxygen().GetVal(cell_id % world.GetWorldX(), cell_id / world.GetWorldX()) == config.INITIAL_OXYGEN_LEVEL()); 
            CHECK(world.GetOxygen().GetNextVal(cell_id % world.GetWorldX(), cell_id / world.GetWorldX()) == Approx(-(config.BASAL_OXYGEN_CONSUMPTION() * config.INITIAL_OXYGEN_LEVEL()/(config.INITIAL_OXYGEN_LEVEL() + config.KM())))); 
        } else {
            CHECK(world.GetOxygen().GetVal(cell_id % world.GetWorldX(), cell_id / world.GetWorldX()) == config.INITIAL_OXYGEN_LEVEL()); 
            CHECK(world.GetOxygen().GetNextVal(cell_id % world.GetWorldX(), cell_id / world.GetWorldX()) == 0); 
        }
    }

    world.GetOxygen().SetVal(example_id % world.GetWorldX(), example_id / world.GetWorldX(),0,0);
    world.RunStep();
    CHECK((!world.IsOccupied(example_id) || world.GetOrg(example_id).hif1alpha == 1));

    world.Reset(config);
    config.OXYGEN_DIFFUSION_COEFFICIENT(.09);
    world.InitConfigs(config);
    CHECK(world.GetOxygen().GetDiffusionCoefficient() == .09);
    world.Run();

}