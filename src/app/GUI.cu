#include "core/Simulation.h"
#include "app/GUI.h"

namespace boltzmann {
    namespace app {
        GUI::GUI(sf::RenderWindow *render_window_, boltzmann::core::Simulation *simulation_)
                : render_window(render_window_), simulation(simulation_) {
            pixels = new sf::Vertex[this->simulation->ydim * this->simulation->xdim];
            vertex_buffer.create(static_cast<size_t>(this->simulation->ydim * this->simulation->xdim));
            cudaMallocManaged(&colors, sizeof(sf::Color) * n_colors);
            cudaMallocManaged(&pixels, sizeof(sf::Vertex) * this->simulation->ydim * this->simulation->xdim);

            for (int y = 0; y < this->simulation->ydim; y++) {
                for (int x = 0; x < this->simulation->xdim; x++) {
                    pixels[y * this->simulation->xdim + x].position = sf::Vector2f{static_cast<float>(x), static_cast<float>(y)};
                }
            }

            sf::Color start_color{0, 0, 0};
            sf::Color end_color{255, 255, 255};
            for (int c = 0; c < n_colors; c++) {
                double percent = (double) c / (double) n_colors;

                auto r = static_cast<uint8_t>((double) start_color.r +
                                              percent * ((double) end_color.r - (double) start_color.r));

                auto g = static_cast<uint8_t>((double) start_color.g +
                                              percent * ((double) end_color.g - (double) start_color.g));

                auto b = static_cast<uint8_t>((double) start_color.b +
                                              percent * ((double) end_color.b - (double) start_color.b));

                colors[c] = sf::Color{r, g, b};
            }
        }

        GUI::~GUI() {
            cudaFree(pixels);
            cudaFree(colors);
        }

        void GUI::paint() {
            #pragma omp parallel for
            for (int y = 0; y < this->simulation->ydim; y++) {
                for (int x = 0; x < this->simulation->xdim; x++) {
                    if (this->simulation->barrier[y][x]) {
                        this->pixels[y * simulation->xdim + x].color = {255, 255, 255};
                    } else {
                        auto colorIndex = std::min(this->n_colors - 1,
                                                   (int) (n_colors * (0.5 + this->simulation->curl[y][x] * contrast * 0.3)));
                        colorIndex = std::max(0, colorIndex);
                        colorIndex = std::min(n_colors - 1, colorIndex);
                        this->pixels[y * simulation->xdim + x].color = this->colors[colorIndex];
                    }
                }
            }
            vertex_buffer.update(pixels);
            render_window->draw(vertex_buffer);
            render_window->display();
        }
    }
}