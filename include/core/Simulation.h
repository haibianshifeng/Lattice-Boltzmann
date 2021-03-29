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

            double** n0_temp;
            double** nN_temp;
            double** nS_temp;
            double** nE_temp;
            double** nW_temp;
            double** nNW_temp;
            double** nNE_temp;
            double** nSW_temp;
            double** nSE_temp;
            double** density_temp;
            double** xvel_temp;
            double** yvel_temp;
            double** speed2_temp;
            double** curl_temp;
            bool** barrier_temp;

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

            void collide() const;

            void stream() const;

            void bounce() const;

            void draw_circle(int x_center, int y_center, int radius) const;

            void compute_curl() const;

            void synchronize() const;

            void debug_information();
        };
    }
}

#endif //BOLTZMANN_SIMULATION_H
