#ifndef _RESOURCE_GRADIENT_H
#define _RESOURCE_GRADIENT_H

#include "base/vector.h"

class ResourceGradient {
    using grid_t = emp::vector<emp::vector<double>>; // It might be possible to make these arrays
    grid_t grid;
    double diffusion_coefficient;
    int x_len;
    int y_len;
    int z_len;

    public:
    ResourceGradient(int x_len_in, int y_len_in=1, int z_len_in=1) :
        x_len(x_len_in), y_len(y_len_in), z_len(z_len_in) {
        for (int y = 0; y < y_len; y++) {
            for (int x = 0; x < x_len; x++) {
                grid[y][x] = 0;
            }
        }
    }

    ResourceGradient(const grid_t & g) {
        x_len = g[0].size();
        y_len = g.size();
        // TODO z_len

        for (int y = 0; y < y_len; y++) {
            for (int x = 0; x < x_len; x++) {
                grid[y][x] = g[y][x];
            }
        }        
    }

    void SetVal(int x, int y, double val) {
        grid[y][x] = val;
    }

    void DecVal(int x, int y, double val) {
        grid[y][x] -= val;
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