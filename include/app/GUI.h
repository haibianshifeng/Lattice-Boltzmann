#ifndef BOLTZMANN_GUI_H
#define BOLTZMANN_GUI_H

#include <SFML/System.hpp>
#include <SFML/Graphics.hpp>
#include <SFML/Audio.hpp>
#include <iostream>
#include <memory>
#include <SFML/OpenGL.hpp>
#include "utils/Exception.h"
#include "core/Simulation.h"
#include "utils/Colors.h"
#include "core/Kernels.h"
#include "utils/Measurement.h"

namespace boltzmann {
    namespace app {
        /**
         * Application gui
         */
        class GUI {
        private:
            /*
             * How many individual colors do we want to have
             */
            int n_colors = 12000;
            sf::Color *colors;

            /*
             * Contrast factor
             */
            double contrast = 400;

            /*
             * Project's specific objects
             */
            boltzmann::core::Simulation *simulation;
            sf::Window *render_window;

            /*
             * OpenGL specific object
             */
            // OpenGL requires a static array of floats to represent coords : [x1][y1][x2][y2]...[xn][yn]
            float **coordinates;
            // OpenGL requires a static array to represent the colors of the dots : [r1][g1][b1][r2][g2][b2]...[rn][gn][bn]
            uint8_t **pixels;
        public:
            /**
             * Constructor
             *
             * @param render_window_ SFML window
             * @param simulation_ simulation object
             */
            GUI(sf::Window *render_window_, core::Simulation *simulation_);

            /**
             * Destructor
             */
            virtual ~GUI();

            /**
             * Visualize current state of the world
             */
            void paint();
        };
    }
}

#endif //BOLTZMANN_GUI_H
