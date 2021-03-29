#include "core/Simulation.h"
#include "app/GUI.h"

namespace boltzmann {
    namespace app {

        __global__
        void update_pixels(uint32_t ydim, uint32_t xdim, sf::Vertex * pixels, bool ** barrier, uint32_t n_colors, double **curl, double contrast, sf::Color * colors) {
            uint32_t y = blockIdx.x;
            uint32_t x = threadIdx.x;

            if(y < ydim && x < xdim) {
                if (barrier[y][x]) {
                    pixels[y * xdim + x].color.r = 255;
                    pixels[y * xdim + x].color.g = 255;
                    pixels[y * xdim + x].color.b = 255;
                } else {
                    auto colorIndex = min(n_colors - 1,
                                          (int) (n_colors *
                                                 (0.5 + curl[y][x] * contrast * 0.3)));
                    colorIndex = max(0, colorIndex);
                    colorIndex = min(n_colors - 1, colorIndex);
                    pixels[y * xdim + x].color = colors[colorIndex];
                }
            }
        }

        GUI::GUI(sf::RenderWindow *render_window_, boltzmann::core::Simulation *simulation_)
                : render_window(render_window_), simulation(simulation_) {
            pixels = new sf::Vertex[this->simulation->ydim * this->simulation->xdim];
            vertex_buffer.create(static_cast<size_t>(this->simulation->ydim * this->simulation->xdim));
            cudaMallocManaged(&colors, sizeof(sf::Color) * n_colors);
            cudaMallocManaged(&pixels, sizeof(sf::Vertex) * this->simulation->ydim * this->simulation->xdim);

            int i = 0;
            for (int y = 0; y < this->simulation->ydim; y++) {
                for (int x = 0; x < this->simulation->xdim; x++) {
                    pixels[i++].position = sf::Vector2f{static_cast<float>(x), static_cast<float>(y)};
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
            update_pixels<<<simulation->ydim, simulation->xdim>>>(simulation->ydim, simulation->xdim, pixels, simulation->barrier, n_colors, simulation->curl, contrast, colors);
            cudaDeviceSynchronize();
            vertex_buffer.update(pixels);
            render_window->draw(vertex_buffer);
            render_window->display();
        }
    }
}