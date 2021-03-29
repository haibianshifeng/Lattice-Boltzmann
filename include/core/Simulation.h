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
#include "core/Kernels.h"

namespace boltzmann {
    namespace core {
        /**
         * Core simulation object
         *
         * Each pixel corresponds a particle
         */
        class Simulation{
        public:
            /*
             * Height and width
             */
            int xdim;
            int ydim;

            /*
             * Main buffer, names and functionalities inspired from the PDF included in this repository
             */
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

            /*
             * Temporary buffer, which has identical content as the main buffer. Since we want to parallelize
             * the streaming and bouncing step (which are difficult to parallelize, we must create those buffers.
             */
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

            /*
             * Speed of particle and omega factor (functionality can be looked at in the PDF included in the repository)
             */
            double v = 0.1;
            double omega = 0.1;

            /**
             * Constructor of class
             *
             * @param width_ Can be any size
             * @param height_ maximal 1024 since CUDA only allows us to use 1024 threads per block.
             */
            Simulation(int width_, int height_, const std::string& barrier_file_name);

            /**
             * Destructor
             */
            virtual ~Simulation();

            /**
             * Initialize begin condition
             */
            void init_fluid() const;

            /**
             * Set buffer at position x and y to 0
             */
            void zero(int x, int y) const;

            /**
             * Draw wall at position x and y
             */
            void draw_barrier(int x, int y) const;

            /**
             * Collision step
             */
            void collide() const;

            /**
             * Stream step
             */
            void stream() const;

            /**
             * Bounce step
             */
            void bounce() const;

            /**
             * Draw a round circle
             *
             * @param x_center x center of circle
             * @param y_center y center of circle
             * @param radius radius of circle
             */
            void draw_circle(int x_center, int y_center, int radius) const;

            /**
             * Compute the curl of the current flow
             */
            void compute_curl() const;

            /**
             * Synchronize temporary buffer with main buffer
             */
            void synchronize() const;

            double getOmega() const;

            void setOmega(double omega);
        };
    }
}

#endif //BOLTZMANN_SIMULATION_H
