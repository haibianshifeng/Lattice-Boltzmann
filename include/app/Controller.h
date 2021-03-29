#ifndef BOLTZMANN_CONTROLLER_H
#define BOLTZMANN_CONTROLLER_H

#include <SFML/System.hpp>
#include <SFML/Graphics.hpp>
#include <SFML/Audio.hpp>
#include <iostream>
#include <memory>
#include "core/Simulation.h"
#include "utils/TimeIt.h"
#include "GUI.h"

namespace boltzmann {
    namespace app {
        /**
         * Application controller
         */
        class Controller {
        private:
            /*
             * SFML's specific objects
             */
            sf::Window *render_window;
            sf::Event event{};

            /*
             * Project's specific objects
             */
            boltzmann::app::GUI *gui;
            boltzmann::core::Simulation * simulation;
        public:

            /**
             * Constructor
             *
             * @param render_window_ SFML main window
             * @param gui_ GUI object
             * @param simulation_ simulation object
             */
            Controller(sf::Window *render_window_, GUI *gui_, core::Simulation *simulation_);

            /**
             * Main loop
             */
            void start(bool recording_mode);
        };
    }
}

#endif //BOLTZMANN_CONTROLLER_H
