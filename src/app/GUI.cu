#include "app/GUI.h"

namespace boltzmann {
    namespace app {
        GUI::GUI(sf::Window *render_window_, boltzmann::core::Simulation *simulation_, bool freaky_colors)
                : render_window(render_window_), simulation(simulation_), colorful(freaky_colors) {
            // Allocate memory for OpenGL coordinates
            cudaMallocManaged(&coordinates, sizeof(float *) * simulation->ydim);

            // Allocate memory for OpenGL colors
            cudaMallocManaged(&pixels, sizeof(uint8_t *) * simulation->ydim);

            // Allocate memory for rainbow colors spectrum
            cudaMallocManaged(&colors, sizeof(sf::Color) * n_colors);

            /*
             * Initialize memory for OpenGL coordinates and colors
             */
            for (int y = 0; y < this->simulation->ydim; y++) {
                cudaMallocManaged(&coordinates[y], sizeof(float) * 2 * simulation->xdim);
                cudaMallocManaged(&pixels[y], sizeof(uint8_t) * 3 * simulation->xdim);

                for (int x = 0; x < this->simulation->xdim; x++) {
                    coordinates[y][2 * x] = (float) x;
                    coordinates[y][2 * x + 1] = (float) y;
                    if (this->simulation->barrier[y][x]) {
                        pixels[y][x * 3] = 50;
                        pixels[y][x * 3 + 1] = 50;
                        pixels[y][x * 3 + 2] = 50;
                    }
                }
            }

            /*
             * Initialize rainbow color for the spectrum
             */
            this->setColorful(freaky_colors);
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

        void GUI::paint(uint32_t mode) {
            /*
             * Rendering the pixels
             * 0: Curl
             * 1: speed
             * 2: x velocity
             * 3: y velocity
             * 4: density
             */
            switch (mode) {
                case 0:
                    boltzmann::core::update_pixels_curl<<<simulation->xdim, simulation->ydim>>>(
                    simulation->ydim,
                            simulation->xdim,
                            pixels,
                            simulation->barrier,
                            n_colors,
                            simulation->curl,
                            contrast,
                            colors);
                    break;
                case 1:
                    boltzmann::core::update_pixels_speed<<<simulation->xdim, simulation->ydim>>>(
                    simulation->ydim,
                            simulation->xdim,
                            pixels,
                            simulation->barrier,
                            n_colors,
                            simulation->speed2,
                            contrast,
                            colors);
                    break;
                case 2:
                    boltzmann::core::update_pixels_xvel<<<simulation->xdim, simulation->ydim>>>(
                    simulation->ydim,
                            simulation->xdim,
                            pixels,
                            simulation->barrier,
                            n_colors,
                            simulation->xvel,
                            contrast,
                            colors);
                    break;
                case 3:
                    boltzmann::core::update_pixels_yvel<<<simulation->xdim, simulation->ydim>>>(
                    simulation->ydim,
                            simulation->xdim,
                            pixels,
                            simulation->barrier,
                            n_colors,
                            simulation->yvel,
                            contrast,
                            colors);
                    break;
                case 4:
                    boltzmann::core::update_pixels_density<<<simulation->xdim, simulation->ydim>>>(
                    simulation->ydim,
                            simulation->xdim,
                            pixels,
                            simulation->barrier,
                            n_colors,
                            simulation->density,
                            contrast,
                            colors);
                    break;
                default:
                    boltzmann::core::update_pixels_curl<<<simulation->xdim, simulation->ydim>>>(
                    simulation->ydim,
                            simulation->xdim,
                            pixels,
                            simulation->barrier,
                            n_colors,
                            simulation->curl,
                            contrast,
                            colors);
            }
            cudaDeviceSynchronize();

            /*
             * Plotting
             */
            glClearColor(0, 0, 0, 0);
            glClear(GL_COLOR_BUFFER_BIT);

            glPushMatrix();
            glEnableClientState(GL_VERTEX_ARRAY);
            glEnableClientState(GL_COLOR_ARRAY);

            for (uint32_t y = 0; y < this->simulation->ydim; y++) {
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

        double GUI::getContrast() const {
            return contrast;
        }

        void GUI::setContrast(double contrast_) {
            GUI::contrast = contrast_;
        }

        bool GUI::isColorful() const {
            return colorful;
        }

        void GUI::setColorful(bool colorful_) {
            this->colorful = colorful_;

            if (this->colorful) {
                for (int c = 0; c < n_colors; c++) {
                    double h = (double) c / n_colors;
                    h += 10 * sin(2 * M_PI * h);
                    colors[c] = boltzmann::utils::HSBtoRGB((float) h, 0.75, 1);
                }
            } else {
                for (int c = 0; c < n_colors; c++) {
                    double h = (2.0 / 3) * (1 - c * 1.0 / n_colors);
                    h += 0.03 * sin(6 * M_PI * h);
                    colors[c] = boltzmann::utils::HSBtoRGB((float) h, (float) 1, (float)1);
                }
            }
        }
    }
}
