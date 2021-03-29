#include <SFML/Graphics.hpp>
#include <app/Controller.h>
#include <SFML/OpenGL.hpp>
#include <app/GUI.h>


int main() {
    constexpr uint32_t width = 1000;
    constexpr uint32_t height = 1000;

    sf::Window window(sf::VideoMode{width, height, 24}, "Boltzmann");
    window.setFramerateLimit(60);

    glViewport(0, 0, width, height); // viewport definition
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0, width, height, 0, -1, 1);
    glPointSize(1);


    boltzmann::core::Simulation simulation(width, height);

    boltzmann::app::GUI gui(&window, &simulation);

    boltzmann::app::Controller controller(&window, &gui, &simulation);

    controller.start();
}