#ifndef _RESOURCE_GRADIENT_H
#define _RESOURCE_GRADIENT_H

#include "base/vector.h"

class ResourceGradient {
    using grid_t = emp::vector<emp::vector<double>>; // It might be possible to make these arrays
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

        curr_grid.resize(y_len);
        next_grid.resize(y_len);

        for (int y = 0; y < y_len; y++) {
            curr_grid[y].resize(x_len);
            next_grid[y].resize(x_len);
            for (int x = 0; x < x_len; x++) {
                curr_grid[y][x] = 0;
                next_grid[y][x] = 0;
            }
        }
    }

    ResourceGradient(const grid_t & g) {
        x_len = g[0].size();
        y_len = g.size();
        // TODO z_len

        curr_grid.resize(y_len);
        next_grid.resize(y_len);

        for (int y = 0; y < y_len; y++) {
            curr_grid[y].resize(x_len);
            next_grid[y].resize(x_len);
            for (int x = 0; x < x_len; x++) {
                curr_grid[y][x] = g[y][x];
                next_grid[y][x] = 0;
            }
        }        
    }

    void SetVal(int x, int y, double val) {
        curr_grid[y][x] = val;
    }

    void SetNextVal(int x, int y, double val) {
        std::cout << "Setting " << x << ", " << y << " to " << val << std::endl;
        next_grid[y][x] = val;
    }

    void DecVal(int x, int y, double val) {
        curr_grid[y][x] -= val;
    }

    void DecNextVal(int x, int y, double val) {
        next_grid[y][x] -= val;
    }

    double GetVal(int x, int y) const {
        return curr_grid[y][x];
    } 

    double GetNextVal(int x, int y) const {
        return next_grid[y][x];
    } 

    void SetDiffusionCoefficient(double coef) {
        diffusion_coefficient = coef;
    }

    void SetToroidal(bool tor) {
        toroidal = tor;
    }

    void Update() {
        std::swap(curr_grid, next_grid);
        for (int y = 0; y < y_len; y++) {
            for (int x = 0; x < x_len; x++) {
                // zero out new next grid
                next_grid[y][x] = 0;

                // Make sure there are no negative numbers in the
                // new curr_grid
                if (curr_grid[y][x] < 0) {
                    curr_grid[y][x] = 0;
                }
            }
        }        
    }

    double GetNeighborOxygen(int x, int y) {
        double total = 0;

        if (toroidal) {
            // Handle left
            if (x <= 0) {
                // Wrap around if toroidal and on left edge
                total += curr_grid[y][x_len-1];
            } else {
                // Otherwise no adjustment is needed
                total += curr_grid[y][x-1];
            }

            // Handle right
            if (x + 1 >= x_len) {
                // Wrap around if toroidal and on right edge
                total += curr_grid[y][0];
            } else {
                total += curr_grid[y][x+1];
            }

            // Handle top
            if (y <= 0) {
                // Wrap around if toroidal and on top edge
                total += curr_grid[y_len-1][x];
            } else {
                // Otherwise no adjustment is needed
                total += curr_grid[y-1][x];
            }

            // Handle bottom
            if (y + 1 >= y_len) {
                // Wrap around if toroidal and on bottom edge
                total += curr_grid[0][x];
            } else {
                total += curr_grid[y+1][x];
            }

        } else {
            // No-flux/Dirichlet

            // Handle left
            if (x <= 0) {
                // Wrap around if toroidal and on left edge
                total += curr_grid[y][x];
            } else {
                // Otherwise no adjustment is needed
                total += curr_grid[y][x-1];
            }

            // Handle right
            if (x + 1 >= x_len) {
                // Wrap around if toroidal and on right edge
                total += curr_grid[y][x];
            } else {
                total += curr_grid[y][x+1];
            }

            // Handle top
            if (y <= 0) {
                // Wrap around if toroidal and on top edge
                total += curr_grid[y][x];
            } else {
                // Otherwise no adjustment is needed
                total += curr_grid[y-1][x];
            }

            // Handle bottom
            if (y + 1 >= y_len) {
                // Wrap around if toroidal and on bottom edge
                total += curr_grid[y][x];
            } else {
                total += curr_grid[y+1][x];
            }
        }

        return total;
    }

    void Diffuse() {

        for (int x = 0; x < x_len; x++) {
            for (int y = 0; y < y_len; y++) {
                next_grid[y][x] += curr_grid[y][x] + 
                        (diffusion_coefficient * 
                        (GetNeighborOxygen(x, y) - 
                        (4.0 * curr_grid[y][x]))); // 4.0 is from central difference approximation
            }
        }
    }


};

#endif