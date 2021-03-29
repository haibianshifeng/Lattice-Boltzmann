#include <SFML/Graphics.hpp>
#include <SFML/OpenGL.hpp>
#include <iostream>
#include "utils/TimeIt.h"

using namespace std;
using namespace sf;

int main() {
    int width = 1000; // window definition
    int height = 1000;
    Window window(sf::VideoMode(width, height, 32), "OpenGL particles");
    window.setFramerateLimit(60);

    glViewport(0, 0, width, height); // viewport definition
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0, width, height, 0, -1, 1);


    auto **coordinates = new float *[height];
    for (int y = 0; y < height; y++) {
        coordinates[y] = new float[2 *width];
        for(int x = 0; x < width; x++) {
            coordinates[y][x * 2] = x;
            coordinates[y][x * 2 + 1] = y;
        }
    }

    auto **colors = new uint8_t *[height];
    for (int y = 0; y < height; y++) {
        colors[y] = new uint8_t[3 * width];
        for(int x = 0; x < width; x++) {
            colors[y][x * 3] = 255;
            colors[y][x * 3 + 1] = 125;
            colors[y][x * 3 + 2] = 125;
        }
    }

    glEnable(GL_POINT_SMOOTH);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glPointSize(1);

    while (window.isOpen()) {
        boltzmann::utils::TimeIt timeIt("Rendering");
        Event event{};
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed)
                window.close();
        }

        glClearColor(0, 0, 0, 0);
        glClear(GL_COLOR_BUFFER_BIT);

        glPushMatrix();
        glEnableClientState(GL_VERTEX_ARRAY);
        glEnableClientState(GL_COLOR_ARRAY);

        for (int y = 0; y < height; y++) {
            glVertexPointer(2, GL_FLOAT, 0, coordinates[y]);
            glColorPointer(3, GL_UNSIGNED_BYTE, 0, colors[y]);
            glDrawArrays(GL_POINTS, 0, width);
        }

        glDisableClientState(GL_VERTEX_ARRAY);
        glDisableClientState(GL_COLOR_ARRAY);

        glPopMatrix();

        glFlush();
        window.display();
        timeIt.end();
    }

    return 0;
}