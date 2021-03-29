#ifndef BOLTZMANN_TIMEIT_H
#define BOLTZMANN_TIMEIT_H

#include <string>
#include <sstream>
#include <chrono>
#include <iostream>
#include <algorithm>

namespace boltzmann {
    namespace utils {
        /**
         * Benchmarking class
         */
        class TimeIt {
        private:
            std::chrono::steady_clock::time_point start_time;
            static uint32_t recursion_level;
            static bool closed_last_line;
        public:
            /**
             * Constructor
             * @param text name of the benchmark
             */
            explicit TimeIt(std::string text, bool verbose);


            /**
             * End the benchmark and print the result
             */
            void end(bool verbose);

            static std::string format_number(uint64_t number);

            static void indent();
        };
    }
}

#endif // BOLTZMANN_TIMEIT_H
