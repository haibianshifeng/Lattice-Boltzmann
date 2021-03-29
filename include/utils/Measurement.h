#ifndef LATTICE_BOLTZMANN_MEASUREMENT_H
#define LATTICE_BOLTZMANN_MEASUREMENT_H
#include <SFML/Graphics.hpp>

namespace boltzmann {
    namespace utils {
        class FPS {
        public:
            FPS() : mFrame(0), mFps(0) {}

            [[nodiscard]] unsigned int getFPS() const { return mFps; }

        private:
            unsigned int mFrame;
            unsigned int mFps;
            sf::Clock mClock;

        public:
            void update() {
                if (mClock.getElapsedTime().asSeconds() >= 1.f) {
                    mFps = mFrame;
                    mFrame = 0;
                    mClock.restart();
                }

                ++mFrame;
            }
        };
    }
}
#endif //LATTICE_BOLTZMANN_MEASUREMENT_H
