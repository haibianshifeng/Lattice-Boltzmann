#include <SFML/Graphics.hpp>
#include <app/Controller.h>
#include <app/GUI.h>


int main() {
    constexpr uint32_t width = 500;
    constexpr uint32_t height = 500;

    sf::RenderWindow window({width, height}, "Boltzmann");

    boltzmann::core::Simulation simulation(width, height);

    boltzmann::app::GUI gui(&window, &simulation);

    boltzmann::app::Controller controller(&window, &gui, &simulation);

    controller.start();
}