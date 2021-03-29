#ifndef BOLTZMANN_GUI_H
#define BOLTZMANN_GUI_H

#include <SFML/System.hpp>
#include <SFML/Graphics.hpp>
#include <SFML/Audio.hpp>
#include <iostream>
#include <memory>
#include "utils/Exception.h"

namespace boltzmann {
    namespace app {
        class GUI {
        private:
            sf::RenderWindow *render_window;
            boltzmann::core::Simulation * simulation;
            sf::Vertex * pixels;
            sf::VertexBuffer vertex_buffer;
            int n_colors = 1200;
            sf::Color * colors;
            double contrast = 400;
        public:
            GUI(sf::RenderWindow *render_window_, core::Simulation *simulation_);

            virtual ~GUI();

            void paint();
        };
    }
}

#endif //BOLTZMANN_GUI_H
