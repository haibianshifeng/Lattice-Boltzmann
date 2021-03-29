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
        class Controller {
        private:
            sf::Window *render_window;
            boltzmann::app::GUI *gui;
            boltzmann::core::Simulation * simulation;
            sf::Event event{};
        public:
            Controller(sf::Window *render_window_, GUI *gui_, core::Simulation *simulation_);

            void start();
        };
    }
}

#endif //BOLTZMANN_CONTROLLER_H
