#include "app/Controller.h"

namespace boltzmann {
    namespace app {
        Controller::Controller(sf::Window *render_window_, GUI *gui_, core::Simulation *simulation_, bool verbose_) : render_window(
                render_window_), gui(gui_), simulation(simulation_) , verbose(verbose_){}

        void Controller::start(bool recording_mode) {
            uint32_t mode = 0;
            bool rendering = false;

            while (this->render_window->isOpen()) {
                while (this->render_window->pollEvent(event)) {
                    switch (this->event.type) {
                        case sf::Event::Closed:
                            this->render_window->close();
                            break;
                        case sf::Event::KeyPressed:
                            if (this->event.key.code == sf::Keyboard::Num0) {
                                mode = 0;
                            } else if (this->event.key.code == sf::Keyboard::Num1) {
                                mode = 1;
                            } else if (this->event.key.code == sf::Keyboard::Num2) {
                                mode = 2;
                            } else if (this->event.key.code == sf::Keyboard::Num3) {
                                mode = 3;
                            } else if (this->event.key.code == sf::Keyboard::Num4) {
                                mode = 4;
                            } else if (this->event.key.code == sf::Keyboard::Add) {
                                this->gui->setContrast(this->gui->getContrast() + 100);
                            } else if (this->event.key.code == sf::Keyboard::Subtract) {
                                this->gui->setContrast(this->gui->getContrast() - 100);
                            } else if(this->event.key.code == sf::Keyboard::Up) {
                                this->simulation->setOmega(this->simulation->getOmega() + 0.001);
                            } else if(this->event.key.code == sf::Keyboard::Down) {
                                this->simulation->setOmega(this->simulation->getOmega() - 0.001);
                            }
                            break;
                        case sf::Event::MouseButtonPressed:
                            if(this->event.mouseButton.button == sf::Mouse::Left) {
                                if (recording_mode) {
                                    if (rendering) {
                                        rendering = false;
                                    } else {
                                        rendering = true;
                                    }
                                }
                            }else if(this->event.mouseButton.button == sf::Mouse::Right){
                                this->gui->setColorful(!this->gui->isColorful());
                            }
                            break;
                        default:
                            break;
                    }
                }
                if(recording_mode && !rendering) {
                    continue;
                } else {

                    if (verbose) {
                        // We simulate four steps at once to make the simulation faster
                        // The number of simulation steps at one time should be even
                        for (int i = 0; i < 4; i++) {
                            boltzmann::utils::TimeIt collision_step("Collision step");
                            this->simulation->collide();
                            cudaDeviceSynchronize();
                            collision_step.end();

                            boltzmann::utils::TimeIt synchronize_step("Synchronization step 1");
                            this->simulation->synchronize();
                            cudaDeviceSynchronize();
                            synchronize_step.end();

                            boltzmann::utils::TimeIt stream_step("Stream step");
                            this->simulation->stream();
                            cudaDeviceSynchronize();
                            stream_step.end();

                            boltzmann::utils::TimeIt synchronize_step1("Synchronization step 2");
                            this->simulation->synchronize();
                            synchronize_step1.end();

                            boltzmann::utils::TimeIt bouncing_step("Bouncing step");
                            this->simulation->bounce();
                            cudaDeviceSynchronize();
                            bouncing_step.end();
                        }

                        // Calculate flow curl at the moment
                        boltzmann::utils::TimeIt curl_step("Curl step");
                        this->simulation->compute_curl();
                        cudaDeviceSynchronize();
                        curl_step.end();

                        // Render world's current state
                        boltzmann::utils::TimeIt rendering_step("Rendering step");
                        this->gui->paint(mode);
                        rendering_step.end();
                    } else {
                        // We simulate four steps at once to make the simulation faster
                        // The number of simulation steps at one time should be even
                        for (int i = 0; i < 4; i++) {
                            this->simulation->collide();
                            cudaDeviceSynchronize();

                            this->simulation->synchronize();
                            cudaDeviceSynchronize();

                            this->simulation->stream();
                            cudaDeviceSynchronize();

                            this->simulation->synchronize();

                            this->simulation->bounce();
                            cudaDeviceSynchronize();
                        }

                        this->simulation->compute_curl();
                        cudaDeviceSynchronize();

                        this->gui->paint(mode);
                    }
                }
            }
        }
    }
}
