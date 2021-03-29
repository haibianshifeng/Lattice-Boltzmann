#include "core/Simulation.h"
#include "app/GUI.h"
#include "utils/Colors.h"

namespace boltzmann {
    namespace app {
        __global__
        void
        update_pixels(uint32_t ydim, uint32_t xdim, sf::Vertex *pixels, bool **barrier, double n_colors, double **curl,
                      double contrast, sf::Color *colors) {
            uint32_t x = blockIdx.x;
            uint32_t y = threadIdx.x;

            if (y < ydim && x < xdim) {
                if (!barrier[y][x]) {
                    auto colorIndex = min(n_colors - 1,
                                          (n_colors *
                                           (0.5f + curl[y][x] * contrast * 0.3f)));
                    colorIndex = max(0.0f, colorIndex);
                    colorIndex = min(n_colors - 1, colorIndex);
                    pixels[y * xdim + x].color = colors[(uint32_t) colorIndex];
                }
            }
        }


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


        }

        GUI::~GUI() {
            cudaFree(pixels);
            cudaFree(colors);
        }

        void GUI::paint() {
            update_pixels<<<simulation->xdim, simulation->ydim>>>(
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
            render_window->display();

        }
    }
}