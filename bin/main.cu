/**
    A lattice-Boltzmann simulation in CUDA

	Credits:
    2011-2013, Daniel V. Schroeder (Weber State University)
	The "wind tunnel" entry/exit conditions are inspired by Graham Pullan's code (http://www.many-core.group.cam.ac.uk/projects/LBdemo.shtml).
    Additional inspiration from Thomas Pohl's applet (http://thomas-pohl.info/work/lba.html).
    Other portions of code are based on Wagner (http://www.ndsu.edu/physics/people/faculty/wagner/lattice_boltzmann_codes/)
    and Gonsalves (http://www.physics.buffalo.edu/phy411-506-2004/index.html;
    code adapted from Succi, http://global.oup.com/academic/product/the-lattice-boltzmann-equation-9780199679249).
	For related materials see:  http://physics.weber.edu/schroeder/fluids
*/

#include <SFML/Graphics.hpp>
#include <app/Controller.h>
#include <SFML/OpenGL.hpp>
#include <app/GUI.h>


int main(int argc, char ** argv) {
    bool recording_mode = false;
    if(argc > 1 && strcmp(argv[1], "record") == 0) {
        recording_mode = true;
    }
    /*
     * Initialize SFML objects for global usage
     */
    constexpr uint32_t width = 1000;
    constexpr uint32_t height = 1000;
    sf::Window window(sf::VideoMode{width, height, 24}, "Boltzmann");
    window.setFramerateLimit(60);

    /*
     * Initialize OpenGL
     */
    glViewport(0, 0, width, height); // viewport definition
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0, width, height, 0, -1, 1);
    glPointSize(1);

    /*
     * Initialize project's specific objects
     */
    boltzmann::core::Simulation simulation(width, height);
    boltzmann::app::GUI gui(&window, &simulation);
    boltzmann::app::Controller controller(&window, &gui, &simulation);

    /*
     * Start main loop
     */
    controller.start(recording_mode);
}