#pragma once

#include "lab3.hpp"

void lab3_start() {
    double time;
    double min_time;
    double max_time;
    double avg_time;

    double res;

    int numbExp = 10;

    min_time = max_time = avg_time = experiment(&res);

    for (int i = 0; i < numbExp; i ++) {
        time = experiment(&res);
        time += res;
        if (max_time < time) max_time = time;
        if (min_time > time) min_time = time;
    }

    std::cout << "execution time: " << avg_time / numbExp << std::endl;
    std::cout << "minimum time: " << min_time << std::endl;
    std::cout << "maximum time: " << max_time << std::endl;
    std::cout << "integral value: " << res << std::endl;
}