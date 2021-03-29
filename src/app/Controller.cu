#include "app/Controller.h"

namespace boltzmann {
    namespace app {
        Controller::Controller(sf::RenderWindow *render_window_, GUI *gui_, core::Simulation *simulation_) : render_window(
                render_window_), gui(gui_), simulation(simulation_) {}

        void Controller::start() {
            while (this->render_window->isOpen()) {
                while (this->render_window->pollEvent(event)) {
                    switch (this->event.type) {
                        case sf::Event::Closed:
                            this->render_window->close();
                            break;
                        default:
                            break;
                    }
                }

                for(int i = 0; i < 10; i++) {
                    boltzmann::utils::TimeIt collision_step("Collision step");
                    this->simulation->collide();
                    collision_step.end();
                    cudaDeviceSynchronize();

                    boltzmann::utils::TimeIt synchronize_step("Synchronization step 1");
                    this->simulation->synchronize();
                    synchronize_step.end();
                    cudaDeviceSynchronize();

                    boltzmann::utils::TimeIt stream_step("Stream step");
                    this->simulation->stream();
                    stream_step.end();
                    cudaDeviceSynchronize();

                    boltzmann::utils::TimeIt synchronize_step1("Synchronization step 2");
                    this->simulation->synchronize();
                    synchronize_step1.end();

                    boltzmann::utils::TimeIt bouncing_step("Bouncing step");
                    this->simulation->bounce();
                    bouncing_step.end();
                }


                boltzmann::utils::TimeIt curl_step("Curl step");
                this->simulation->compute_curl();
                curl_step.end();

                cudaDeviceSynchronize();

                boltzmann::utils::TimeIt rendering_step("Rendering step");
                this->gui->paint();
                rendering_step.end();
            }
        }
    }
}
