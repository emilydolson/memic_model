#ifndef _RESOURCE_GRADIENT_H
#define _RESOURCE_GRADIENT_H

#include "base/array.h"

template <int X_LEN, int Y_LEN = 0, int Z_LEN = 0>
class ResourceGradient {
    using grid_t = emp::array<emp::array<double, X_LEN>, Y_LEN>;
    grid_t grid;
    double diffusion_coefficient;

    ResourceGradient() {
        for (int y = 0; y < Y_LEN; y++) {
            for (int x = 0; x < X_LEN; x++) {
                grid[y][x] = 0;
            }
        }
    }

    ResourceGradient(grid_t g) {
        for (int y = 0; y < Y_LEN; y++) {
            for (int x = 0; x < X_LEN; x++) {
                grid[y][x] = g[y][x];
            }
        }        
    }

    void SetVal(int x, int y, double val) {
        grid[y][x] = val;
    }

    double GetVal(int x, int y) {
        return grid[y][x];
    }

    void Diffuse() {
        // TODO
        return;
    }


};

#endif