#include "app/GUI.h"

namespace boltzmann {
    namespace app {
        GUI::GUI(sf::RenderWindow *render_window_, boltzmann::core::Simulation *simulation_)
                : render_window(render_window_), simulation(simulation_) {
            pixels = new sf::Vertex[this->simulation->ydim * this->simulation->xdim];
            vertex_buffer.create(this->simulation->ydim * this->simulation->xdim);
            cudaMallocManaged(&colors, sizeof(sf::Color) * n_colors);
            cudaMallocManaged(&pixels, sizeof(sf::Vertex) * this->simulation->ydim * this->simulation->xdim);

            for (int y = 0; y < this->simulation->ydim; y++) {
                for (int x = 0; x < this->simulation->xdim; x++) {
                    pixels[y * this->simulation->xdim + x].position = sf::Vector2f{static_cast<float>(x),
                                                                                   static_cast<float>(y)};
                    if(this->simulation->barrier[y][x]) {
                        pixels[y * this->simulation->xdim + x].color = sf::Color{125, 125, 125};
                    }
                }
            }

            for (int c = 0; c < n_colors; c++) {
                double h = (double)c / n_colors;
                h += 3 * sin(4*M_PI*h);
                colors[c] = HSBtoRGB((float) h, 0.75, 1);
            }

            if(!font.loadFromFile("../../data/arial.ttf")) {
                THROW_EXCEPTION("Can not find 'arial.ttf'. Exit now!")
            }
        }

        GUI::~GUI() {
            cudaFree(pixels);
            cudaFree(colors);
        }

        void GUI::paint() {
            boltzmann::core::update_pixels<<<simulation->xdim, simulation->ydim>>>(
                    simulation->ydim,
                    simulation->xdim,
                    pixels,
                    simulation->barrier,
                    n_colors,
                    simulation->curl,
                    contrast,
                    colors);
            cudaDeviceSynchronize();
            vertex_buffer.update(pixels);
            render_window->draw(vertex_buffer);
            this->draw_fps();
            render_window->display();

        }

        void GUI::draw_fps() {
            std::stringstream ss;
            fps_measurement.update();
            ss << "FPS: " << fps_measurement.getFPS();
            sf::Text fpsText{ss.str(), font, 10};
            fpsText.setPosition(simulation->xdim - 50, 10);
            fpsText.setFillColor(sf::Color::Black);
            render_window->draw(fpsText);
        }
    }
}