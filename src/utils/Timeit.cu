#include <utility>

#include "utils/TimeIt.h"

namespace boltzmann {
    namespace utils {
        uint32_t TimeIt::recursion_level;

        bool TimeIt::closed_last_line;

        std::string pad_to(const std::string &str, const size_t num, const char paddingChar = ' ') {
            std::string res = str;
            if (num > res.size())
                res.append(num - res.size(), paddingChar);
            return res;
        }

        TimeIt::TimeIt(std::string text, bool verbose){
            if(verbose ) {
                this->start_time = std::chrono::steady_clock::now();

                if (!TimeIt::closed_last_line)
                    std::cout << std::endl;
                TimeIt::indent();
                std::cout << pad_to(text, 75 - TimeIt::recursion_level * 2) << "  ";

                TimeIt::recursion_level++;
                TimeIt::closed_last_line = false;
            }
        }

        void TimeIt::end(bool verbose) {
            if(verbose) {
                auto end_time = std::chrono::steady_clock::now();

                TimeIt::recursion_level--;

                auto diff = end_time - this->start_time;

                if (TimeIt::closed_last_line) {
                    TimeIt::indent();

                    std::cout << pad_to("- ", 75 - TimeIt::recursion_level * 2) << "  ";
                }
                std::cout << "Done after "
                          << TimeIt::format_number(
                                  static_cast<uint64_t>(std::chrono::duration_cast<std::chrono::microseconds>(
                                          diff).count()))
                          << "[us]" << std::endl;

                TimeIt::closed_last_line = true;
            }
        }

        void TimeIt::indent() {
            for (uint32_t i = 0; i < TimeIt::recursion_level; i++) {
                std::cout << "  ";
            }
        }

        std::string TimeIt::format_number(uint64_t number) {
            auto src = std::to_string(number);
            auto dest = std::string();

            auto count = 3;
            for (auto i = src.crbegin(); i != src.crend(); ++i) {
                if (count == 0) {
                    dest.push_back(',');
                    count = 3;
                }
                if (count--) {
                    dest.push_back(*i);
                }
            }
            std::reverse(dest.begin(), dest.end());
            return dest;
        }
    }
}