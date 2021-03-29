#ifndef BOLTZMANN_SIMULATION_H
#define BOLTZMANN_SIMULATION_H

#include <SFML/System.hpp>
#include <SFML/Graphics.hpp>
#include <SFML/Audio.hpp>
#include <memory>
#include <omp.h>
#include <cmath>
#include <cstring>

#include "utils/Exception.h"

namespace boltzmann {
    namespace core {
        class Simulation{
        public:
            int xdim;
            int ydim;
            double** n0;
            double** nN;
            double** nS;
            double** nE;
            double** nW;
            double** nNW;
            double** nNE;
            double** nSW;
            double** nSE;

            double** density;		
            double** xvel;
            double** yvel;			
            double** speed2;
            double** curl;

            bool** barrier;

            const double four9ths = 4.0 / 9;
            const double one9th = 1.0 / 9;
            const double one36th = 1.0 / 36;
            double v = 0.1;
            double omega = 1.78;

            void init_fluid() const;

            void zeroSite(int x, int y) const;

            Simulation(int width_, int height_);

            virtual ~Simulation();

            void draw_barrier(int x, int y) const;

            void collide();

            void stream() const;

            void bounce() const;

            void draw_circle(int x_center, int y_center, int radius) const;

            void draw_square(int x1, int y1, int x2, int y2) const;

            void compute_curl() const;
        };
    }
}

#endif //BOLTZMANN_SIMULATION_H
