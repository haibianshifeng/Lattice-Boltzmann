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

                for(int i = 0; i < 20; i++) {
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


                boltzmann::utils::TimeIt curl_step("Curl step");
                this->simulation->compute_curl();
                cudaDeviceSynchronize();
                curl_step.end();


                boltzmann::utils::TimeIt rendering_step("Rendering step");
                this->gui->paint();
                cudaDeviceSynchronize();
                rendering_step.end();

            }
        }
    }
}
