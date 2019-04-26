#ifndef _RESOURCE_GRADIENT_H
#define _RESOURCE_GRADIENT_H

#include "base/vector.h"

class ResourceGradient {
    using grid_t = emp::vector<emp::vector<emp::vector<double>>>; // It might be possible to make these arrays
    grid_t curr_grid;
    grid_t next_grid;
    double diffusion_coefficient;
    int x_len;
    int y_len;
    int z_len;
    bool toroidal;

    public:
    ResourceGradient(int x_len_in, int y_len_in=1, int z_len_in=1) :
        diffusion_coefficient(0),
        x_len(x_len_in), y_len(y_len_in), z_len(z_len_in),
        toroidal(false) {

        curr_grid.resize(z_len);
        next_grid.resize(z_len);

        for (int z = 0; z < z_len; z++) {
            curr_grid[z].resize(y_len);
            next_grid[z].resize(y_len);
            for (int y = 0; y < y_len; y++) {
                curr_grid[z][y].resize(x_len);
                next_grid[z][y].resize(x_len);
                for (int x = 0; x < x_len; x++) {
                    curr_grid[z][y][x] = 0;
                    next_grid[z][y][x] = 0;
                }
            }
        }
    }

    ResourceGradient(const grid_t & g) {
        x_len = g[0][0].size();        
        y_len = g[0].size();
        z_len = g.size();        

        curr_grid.resize(z_len);
        next_grid.resize(z_len);

        for (int z = 0; z < z_len; z++) {
            curr_grid[z].resize(y_len);
            next_grid[z].resize(y_len);
            for (int y = 0; y < y_len; y++) {
                curr_grid[z][y].resize(x_len);
                next_grid[z][y].resize(x_len);
                for (int x = 0; x < x_len; x++) {
                    curr_grid[z][y][x] = g[z][y][x];
                    next_grid[z][y][x] = 0;
                }
            }
        }
 
    }

    void SetVal(int x, int y, int z, double val) {
        curr_grid[z][y][x] = val;
    }

    void SetNextVal(int x, int y, int z, double val) {
        std::cout << "Setting " << x << ", " << y << " to " << val << std::endl;
        next_grid[z][y][x] = val;
    }

    void DecVal(int x, int y, int z, double val) {
        curr_grid[z][y][x] -= val;
    }

    void DecNextVal(int x, int y, int z, double val) {
        next_grid[z][y][x] -= val;
    }

    double GetVal(int x, int y, int z=0) const {
        return curr_grid[z][y][x];
    } 

    double GetNextVal(int x, int y, int z = 0) const {
        return next_grid[z][y][x];
    } 

    void SetDiffusionCoefficient(double coef) {
        diffusion_coefficient = coef;
    }

    void SetToroidal(bool tor) {
        toroidal = tor;
    }

    void Update() {
        std::swap(curr_grid, next_grid);
        for (int z = 0; z < z_len; z++) {
            for (int y = 0; y < y_len; y++) {
                for (int x = 0; x < x_len; x++) {
                    // zero out new next grid
                    next_grid[z][y][x] = 0;

                    // Make sure there are no negative numbers in the
                    // new curr_grid
                    if (curr_grid[z][y][x] < 0) {
                        curr_grid[z][y][x] = 0;
                    }
                }
            }
        }        
    }

    double GetNeighborOxygen(int x, int y, int z) {
        double total = 0;

        if (toroidal) {
            // Handle left
            if (x <= 0) {
                // Wrap around if toroidal and on left edge
                total += curr_grid[z][y][x_len-1];
            } else {
                // Otherwise no adjustment is needed
                total += curr_grid[z][y][x-1];
            }

            // Handle right
            if (x + 1 >= x_len) {
                // Wrap around if toroidal and on right edge
                total += curr_grid[z][y][0];
            } else {
                total += curr_grid[z][y][x+1];
            }

            // Handle top
            if (y <= 0) {
                // Wrap around if toroidal and on top edge
                total += curr_grid[z][y_len-1][x];
            } else {
                // Otherwise no adjustment is needed
                total += curr_grid[z][y-1][x];
            }

            // Handle bottom
            if (y + 1 >= y_len) {
                // Wrap around if toroidal and on bottom edge
                total += curr_grid[z][0][x];
            } else {
                total += curr_grid[z][y+1][x];
            }

            // Handle below
            if (z <= 0) {
                // Wrap around if toroidal and on top edge
                total += curr_grid[z_len - 1][y][x];
            } else {
                // Otherwise no adjustment is needed
                total += curr_grid[z-1][y][x];
            }

            // Handle above
            if (z + 1 >= z_len) {
                // Wrap around if toroidal and on bottom edge
                total += curr_grid[0][y][x];
            } else {
                total += curr_grid[z+1][y][x];
            }


        } else {
            // No-flux/Dirichlet

            // Handle left
            if (x <= 0) {
                // Do the Drichelet thing
                total += curr_grid[z][y][x];
            } else {
                // Otherwise no adjustment is needed
                total += curr_grid[z][y][x-1];
            }

            // Handle right
            if (x + 1 >= x_len) {
                // Do the Drichelet thing
                total += curr_grid[z][y][x];
            } else {
                total += curr_grid[z][y][x+1];
            }

            // Handle top
            if (y <= 0) {
                // Do the Drichelet thing
                total += curr_grid[z][y][x];
            } else {
                // Otherwise no adjustment is needed
                total += curr_grid[z][y-1][x];
            }

            // Handle bottom
            if (y + 1 >= y_len) {
                // Do the Drichelet thing
                total += curr_grid[z][y][x];
            } else {
                total += curr_grid[z][y+1][x];
            }

            // Handle below
            if (z <= 0) {
                // Do the Drichelet thing
                total += curr_grid[z][y][x];
            } else {
                // Otherwise no adjustment is needed
                total += curr_grid[z-1][y][x];
            }

            // Handle above
            if (z + 1 >= z_len) {
                // Do the Drichelet thing
                total += curr_grid[z][y][x];
            } else {
                total += curr_grid[z+1][y][x];
            }

        }

        return total;
    }

    void Diffuse() {

        for (int z = 0; z < z_len; z++) {
            for (int x = 0; x < x_len; x++) {
                for (int y = 0; y < y_len; y++) {
                    next_grid[z][y][x] += curr_grid[z][y][x] + 
                            (diffusion_coefficient * 
                            (GetNeighborOxygen(x, y, z) - 
                            (6.0 * curr_grid[z][y][x]))); // 6.0 is from central difference approximation
                }
            }
        }
    }


};

#endif