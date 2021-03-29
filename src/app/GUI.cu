#include "app/GUI.h"

namespace boltzmann {
    namespace app {
        GUI::GUI(sf::Window *render_window_, boltzmann::core::Simulation *simulation_)
                : render_window(render_window_), simulation(simulation_) {
            cudaMallocManaged(&coordinates, sizeof(float *) * simulation->ydim);
            cudaMallocManaged(&pixels, sizeof(uint8_t *) * simulation->ydim);
            cudaMallocManaged(&colors, sizeof(sf::Color) * n_colors);

            for (int y = 0; y < this->simulation->ydim; y++) {
                cudaMallocManaged(&coordinates[y], sizeof(float) * 2 * simulation->xdim);
                cudaMallocManaged(&pixels[y], sizeof(uint8_t) * 3 * simulation->xdim);

                for (int x = 0; x < this->simulation->xdim; x++) {
                    coordinates[y][2 * x] = (float) x;
                    coordinates[y][2 * x + 1] = (float) y;
                    if (this->simulation->barrier[y][x]) {
                        pixels[y][x * 3] = 125;
                        pixels[y][x * 3 + 1] = 125;
                        pixels[y][x * 3 + 2] = 125;
                    }
                }
            }

            for (int c = 0; c < n_colors; c++) {
                double h = (double) c / n_colors;
                h += 3 * sin(4 * M_PI * h);
                colors[c] = HSBtoRGB((float) h, 0.75, 1);
            }

            if (!font.loadFromFile("../../data/arial.ttf")) {
                THROW_EXCEPTION("Can not find 'arial.ttf'. Exit now!")
            }
        }

        GUI::~GUI() {
            for (int y = 0; y < this->simulation->ydim; y++) {
                cudaFree(&coordinates[y]);
                cudaFree(&pixels[y]);
            }
            cudaFree(pixels);
            cudaFree(colors);
            cudaFree(coordinates);
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

            glClearColor(0, 0, 0, 0);
            glClear(GL_COLOR_BUFFER_BIT);

            glPushMatrix();
            glEnableClientState(GL_VERTEX_ARRAY);
            glEnableClientState(GL_COLOR_ARRAY);

            for(uint32_t y = 0; y < this->simulation->ydim; y++) {
                glVertexPointer(2, GL_FLOAT, 0, coordinates[y]);
                glColorPointer(3, GL_UNSIGNED_BYTE, 0, pixels[y]);
                glDrawArrays(GL_POINTS, 0, simulation->xdim);
            }

            glDisableClientState(GL_VERTEX_ARRAY);
            glDisableClientState(GL_COLOR_ARRAY);

            glPopMatrix();
            glFlush();

            render_window->display();
        }
    }
}