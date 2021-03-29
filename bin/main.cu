/**
    A lattice-Boltzmann simulation in CUDA

	Full Credit goes to the creator of the exercises:
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
#include <app/CLI11.hpp>

typedef struct cli_paramters {
    bool recording;
    bool freaky_colors;
    bool verbose;
    uint32_t width;
    uint32_t height;
    std::string barrier_file_name;
}CliParameters;

int main(int argc, char ** argv) {
    /*
     * Parsing command line arguments
     */
    CliParameters  cli_parameters{
        .recording = false,
        .freaky_colors = false,
        .verbose = false,
        .width = 1000,
        .height = 1000
    };
    CLI::App app{"Lattice Boltzmann Simulation"};
    app.add_flag("-r,--recording",
                 cli_parameters.recording,
                 "Record mode on. Default false. If this flag is true, left mouse click on the window is needed to start the simulation.");

    app.add_flag("-f,--freaky",
                 cli_parameters.freaky_colors,
                 "Freaky colors on. Default false. If this flag is true, non-traditional colors will be used, else traditional colors.");

    app.add_flag("-v,--verbose", cli_parameters.verbose, "Verbosity for benchmarking. Default false.");

    app.add_option("-b,--barrier", cli_parameters.barrier_file_name, "Path to png/jpeg/jpg images to import self-made barrier mask file. Darker areas of the image (average RGB less than 100) will be detected as barrier.");

    app.add_option("-x,--width", cli_parameters.width, "Width of the application. (Default 1000).");

    app.add_option("-y,--height", cli_parameters.height, "Height of the application. (Default 1000). Should be at most 1000.");

    CLI11_PARSE(app, argc, argv)

    sf::Window window(sf::VideoMode{cli_parameters.width, cli_parameters.height, 24}, "Boltzmann");
    window.setFramerateLimit(60);

    /*
     * Initialize OpenGL
     */
    glViewport(0, 0, cli_parameters.width, cli_parameters.height); // viewport definition
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0, cli_parameters.width, cli_parameters.height, 0, -1, 1);
    glPointSize(1);

    /*
     * Initialize project's specific objects
     */
    boltzmann::core::Simulation simulation(cli_parameters.width, cli_parameters.height, cli_parameters.barrier_file_name);
    boltzmann::app::GUI gui(&window, &simulation, cli_parameters.freaky_colors);
    boltzmann::app::Controller controller(&window, &gui, &simulation, cli_parameters.verbose);

    /*
     * Start main loop
     */
    controller.start(cli_parameters.recording);
}