#ifndef BOLTZMANN_GUI_H
#define BOLTZMANN_GUI_H

#include <SFML/System.hpp>
#include <SFML/Graphics.hpp>
#include <SFML/Audio.hpp>
#include <iostream>
#include <memory>
#include "utils/Exception.h"
#include "core/Simulation.h"
#include "utils/Colors.h"
#include "core/Kernels.h"
#include "utils/Measurement.h"

namespace boltzmann {
    namespace app {
        class GUI {
        private:
            sf::RenderWindow *render_window;
            boltzmann::core::Simulation * simulation;
            sf::Vertex * pixels;
            sf::VertexBuffer vertex_buffer;
            int n_colors = 12000;
            sf::Color * colors;
            double contrast = 400;
            boltzmann::utils::FPS fps_measurement{};
            sf::Font font;
            void draw_fps();
        public:
            GUI(sf::RenderWindow *render_window_, core::Simulation *simulation_);

            virtual ~GUI();

            void paint();
        };
    }
}

#endif //BOLTZMANN_GUI_H
