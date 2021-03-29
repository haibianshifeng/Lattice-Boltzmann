#ifndef LATTICE_BOLTZMANN_KERNELS_H
#define LATTICE_BOLTZMANN_KERNELS_H

#include <SFML/Graphics.hpp>

namespace boltzmann {
    namespace core {
        /*
         * Vectors' probabilities constants
         */
        const double four9ths = 4.0 / 9.0;
        const double one9th = 1.0 / 9.0;
        const double one36th = 1.0 / 36.0;

        /**
         * Collision step
         */
        __global__
        void collide(uint32_t xdim, uint32_t ydim,
                     bool **barrier,
                     double **n0,
                     double **nN,
                     double **nS,
                     double **nE,
                     double **nW,
                     double **nNW,
                     double **nNE,
                     double **nSW,
                     double **nSE,
                     double **density,
                     double **xvel,
                     double **yvel,
                     double **speed2,
                     double **n0_temp,
                     double **nN_temp,
                     double **nS_temp,
                     double **nE_temp,
                     double **nW_temp,
                     double **nNW_temp,
                     double **nNE_temp,
                     double **nSW_temp,
                     double **nSE_temp,
                     double **density_temp,
                     double **xvel_temp,
                     double **yvel_temp,
                     double **speed2_temp,
                     double omega);

        /**
         * Computing flow's curl
         */
        __global__
        void compute_curl(uint32_t xdim, uint32_t ydim, double **curl, double **yvel, double **xvel);

        /**
         * Streaming step
         */
        __global__
        void stream(uint32_t xdim,
                    uint32_t ydim,
                    bool **barrier,
                    double **n0,
                    double **nN,
                    double **nS,
                    double **nE,
                    double **nW,
                    double **nNW,
                    double **nNE,
                    double **nSW,
                    double **nSE,
                    double **density,
                    double **xvel,
                    double **yvel,
                    double **speed2,
                    double **n0_temp,
                    double **nN_temp,
                    double **nS_temp,
                    double **nE_temp,
                    double **nW_temp,
                    double **nNW_temp,
                    double **nNE_temp,
                    double **nSW_temp,
                    double **nSE_temp,
                    double **density_temp,
                    double **xvel_temp,
                    double **yvel_temp,
                    double **speed2_temp,
                    double omega,
                    double v);

        /**
         * Bouncing step
         */
        __global__
        void bounce(uint32_t xdim,
                    uint32_t ydim,
                    bool **barrier,
                    double **n0,
                    double **nN,
                    double **nS,
                    double **nE,
                    double **nW,
                    double **nNW,
                    double **nNE,
                    double **nSW,
                    double **nSE,
                    double **density,
                    double **xvel,
                    double **yvel,
                    double **speed2,
                    double **n0_temp,
                    double **nN_temp,
                    double **nS_temp,
                    double **nE_temp,
                    double **nW_temp,
                    double **nNW_temp,
                    double **nNE_temp,
                    double **nSW_temp,
                    double **nSE_temp,
                    double **density_temp,
                    double **xvel_temp,
                    double **yvel_temp,
                    double **speed2_temp,
                    double omega,
                    double v);

        /**
         * Synchronize main buffer with temporary buffer on graphic card
         */
        __global__
        void synchronize(uint32_t xdim,
                         uint32_t ydim,
                         bool **barrier,
                         double **n0,
                         double **nN,
                         double **nS,
                         double **nE,
                         double **nW,
                         double **nNW,
                         double **nNE,
                         double **nSW,
                         double **nSE,
                         double **density,
                         double **xvel,
                         double **yvel,
                         double **speed2,
                         double **n0_temp,
                         double **nN_temp,
                         double **nS_temp,
                         double **nE_temp,
                         double **nW_temp,
                         double **nNW_temp,
                         double **nNE_temp,
                         double **nSW_temp,
                         double **nSE_temp,
                         double **density_temp,
                         double **xvel_temp,
                         double **yvel_temp,
                         double **speed2_temp,
                         double omega,
                         double v);

        /**
         * Generate visual graphic from current world's state
         */
        __global__
        void
        update_pixels(uint32_t ydim, uint32_t xdim, uint8_t **pixels, bool **barrier, double n_colors, double **curl,
                      double contrast, sf::Color *colors);
    }
}

#endif //LATTICE_BOLTZMANN_KERNELS_H
