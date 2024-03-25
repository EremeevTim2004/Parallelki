#pragma once

#include "lab3.hpp"

double rand_rang(double a, double b) {
    if (a > b) return range_rang(b, a);
    else return a + (b - a) * rand() / RAND_MAX;
}

void integral(const double a1, const double b1, const double a2, const double b2, const double h1, const double h2, double* res) {
    int i = 0;
    int j = 0;

    int n1 = (int)((b1 - a1) / h1);
    int n2 = (int)((b2 - a2) / h2);

    double sum = 0.0;

    double x;
    double y;

#pragma parallel for reduction(+: sum)
    for (i = 1; i <= n1; i ++) {
        x = rand_rang(a1, b1);
        y = rand_rang(a2, b2);
        sum += ((pow(e, sin(PI * x) * cos(PI * y)) + 1) / ((b1 - a1) * (b2 - a2)));
    }
    *res = (b1 - a1) * (b2 - a2) * sum / n;
}

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