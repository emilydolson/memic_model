#define CATCH_CONFIG_MAIN  // This tells Catch to provide a main() - only do this in one cpp file
#include "catch.hpp"
#include "../source/ResourceGradient.h"

int y_len = 50;
int x_len = 100;


TEST_CASE("Test constructor", "[oxygen_gradient]") {
    ResourceGradient r(x_len, y_len);
    
    // Make sure r got constructed correctly
    for (int y = 0; y < y_len; y++) {
        for (int x = 0; x < x_len; x++) {
            CHECK(r.GetVal(x, y) == 0);
        }
    }

    // Test constructing from grid
    emp::vector<emp::vector<double>> grid;
    int grid_y = 5;
    int grid_x = 3;
    grid.resize(grid_y);
    for (int y = 0; y < grid_y; y++) {
        grid[y].resize(grid_x);
        for (int x = 0; x < grid_x; x++) {
            grid[y][x] = x*x+y;
        }
    }

    ResourceGradient r2(grid);
    for (int y = 0; y < grid_y; y++) {
        for (int x = 0; x < grid_x; x++) {
            CHECK(r2.GetVal(x, y) == x*x+y);
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
            CHECK(r.GetVal(x, y) == 0);
            CHECK(r.GetNeighborOxygen(x,y) == 0);
        }
    }
    
    // Add some oxygen (and make sure that it works)
    r.SetVal(0,0,1);
    CHECK(r.GetVal(0,0) == 1);

    // Test GetNeighborOxygen
    // All the oxygen is in 0,0, not its neighbors,
    // But Dirichlet boundaries cause us to assume the cell
    // has two neighbors off-screen with equivalent oxygen
    CHECK(r.GetNeighborOxygen(0,0) == 2);
    CHECK(r.GetNeighborOxygen(1,0) == 1);
    CHECK(r.GetNeighborOxygen(0,1) == 1);
    CHECK(r.GetNeighborOxygen(1,1) == 0);

    // Not toroidal
    CHECK(r.GetNeighborOxygen(99,49) == 0);
    CHECK(r.GetNeighborOxygen(99,0) == 0);
    CHECK(r.GetNeighborOxygen(0,49) == 0);

    // Now do actual diffusion
    r.Diffuse();
    CHECK(r.GetNextVal(0,0) == Approx(0.98));
    CHECK(r.GetNextVal(0,1) == Approx(0.01));
    CHECK(r.GetNextVal(1,0) == Approx(0.01));
    CHECK(r.GetNextVal(1,1) == 0);

    // Not toroidal
    CHECK(r.GetNextVal(99,49) == 0);
    CHECK(r.GetNextVal(99,0) == 0);
    CHECK(r.GetNextVal(0,49) == 0);

}

TEST_CASE("Test turning off diffusion", "[oxygen_gradient]") {
    ResourceGradient r(x_len, y_len);

    // If there is no oxygen, there should still be no oxygen
    // after diffusion happens
    r.SetDiffusionCoefficient(0);
    r.Diffuse();
    for (int y = 0; y < y_len; y++) {
        for (int x = 0; x < x_len; x++) {
            CHECK(r.GetVal(x, y) == 0);
            CHECK(r.GetNeighborOxygen(x,y) == 0);
        }
    }
    
    // Add some oxygen (and make sure that it works)
    r.SetVal(0,0,1);
    CHECK(r.GetVal(0,0) == 1);

    // Test GetNeighborOxygen
    // All the oxygen is in 0,0, not its neighbors
    CHECK(r.GetNeighborOxygen(0,0) == 2);
    CHECK(r.GetNeighborOxygen(1,0) == 1);
    CHECK(r.GetNeighborOxygen(0,1) == 1);
    CHECK(r.GetNeighborOxygen(1,1) == 0);

    // Not toroidal
    CHECK(r.GetNeighborOxygen(99,49) == 0);
    CHECK(r.GetNeighborOxygen(99,0) == 0);
    CHECK(r.GetNeighborOxygen(0,49) == 0);

    // Now do actual diffusion
    r.Diffuse();
    CHECK(r.GetNextVal(0,0) == Approx(1));
    CHECK(r.GetNextVal(0,1) == Approx(0));
    CHECK(r.GetNextVal(1,0) == Approx(0));
    CHECK(r.GetNextVal(1,1) == 0);

    // Not toroidal
    CHECK(r.GetNextVal(99,49) == 0);
    CHECK(r.GetNextVal(99,0) == 0);
    CHECK(r.GetNextVal(0,49) == 0);

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
            CHECK(r.GetVal(x, y) == 0);
            CHECK(r.GetNeighborOxygen(x,y) == 0);
        }
    }
    
    // Add some oxygen (and make sure that it works)
    r.SetVal(0,0,1);
    CHECK(r.GetVal(0,0) == 1);

    // Test GetNeighborOxygen
    // All the oxygen is in 0,0, not its neighbors
    CHECK(r.GetNeighborOxygen(0,0) == 0);
    CHECK(r.GetNeighborOxygen(1,0) == 1);
    CHECK(r.GetNeighborOxygen(0,1) == 1);
    CHECK(r.GetNeighborOxygen(1,1) == 0);

    // Toroidal
    CHECK(r.GetNeighborOxygen(99,49) == 0);
    CHECK(r.GetNeighborOxygen(99,0) == 1);
    CHECK(r.GetNeighborOxygen(0,49) == 1);

    // Now do actual diffusion
    r.Diffuse();
    CHECK(r.GetNextVal(0,0) == Approx(0.96));
    CHECK(r.GetNextVal(0,1) == Approx(0.01));
    CHECK(r.GetNextVal(1,0) == Approx(0.01));
    CHECK(r.GetNextVal(1,1) == 0);

    // Toroidal
    CHECK(r.GetNextVal(99,49) == 0);
    CHECK(r.GetNextVal(99,0) == 0.01);
    CHECK(r.GetNextVal(0,49) == 0.01);

}

TEST_CASE("Test updating gradient", "[oxygen_gradient]") {
    ResourceGradient r(x_len, y_len);

    r.SetVal(10,10,5);
    r.SetVal(1,1,5);
    CHECK(r.GetVal(10,10) == 5);
    CHECK(r.GetVal(1,1) == 5);

    r.DecVal(10,10, 1);
    r.SetNextVal(8,6,3);
    r.SetNextVal(10,10,20);
    r.DecNextVal(10,10,2);
    CHECK(r.GetVal(10,10) == 4);
    CHECK(r.GetNextVal(8,6) == 3);
    CHECK(r.GetNextVal(10,10) == 18);

    r.Update();
    CHECK(r.GetVal(8,6) == 3);
    CHECK(r.GetVal(10,10) == 18);
    CHECK(r.GetVal(1,1) == 0);
}