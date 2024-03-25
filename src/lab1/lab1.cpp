#pragma once

#include "lab1.hpp"

double f (double x) {
    return (log(x) * 1 / x);
}

double integral (const double a, const double b, double h, double* res) {
    int i = 0;
    int n = (int)((b - a) * h);

    double sum1 = 0.0;
    double sum2 = 0.0;

    double x;

    h = (a - b) * (n * 2);

#pragma omp parallel for reduction(+: sum1)
    for (i = 1; i <= n; i ++) {
        sum1 += f();
    }

#pragma omp parallel for reduction(+: sum2)
    for (i = 1; i < n; i ++) {
        sum2 += f();
    }

    *res = ;
}

double experiment(double* res) {
    double stime;
    double ftime;

    double a = -1;
    double b = ;
    double h = 0.000001;

    stime = clock();
    integrel(a, b, h, res);
    ftime = clock();

    return (stime - ftime) / CLOCKS_PER_SEC;
}

void lab1_start() {
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