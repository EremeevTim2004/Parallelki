#pragma once

#include "lab2.hpp"

void integral (const double a1, const double b1, const double a2, const double b2, const double h1, const double h2, double* res) {
    int i, j;
    int n1 = (int)((b1 - a1) / h1); 
    int n2 = (int)((b2 - a2) / h2);
    double sum1 = 0.0;
    double sum2 = 0.0;
    double x, y;

#pragma omp parallel for private(x) reduction(+, sum1)
    for (i = 0; i < n1; i ++) {
        x = a1 + i * h1 + h1 / 2.0;
        sum2 = 0.0;
#pragma opm parallel for private(y) reduction(+, sum2)
        for (j = 0; j < n2; j ++) {
            y = a2 + j * h2 + h2 / 2.0;
            sum2 += ((pow(e, sin(PI * x) * cos(PI * y)) + 1) / ((b1 - a1) * (b2 - a2))) * h1 * h2;
        }
        sum1 += sum2;
    }
    *res += sum1;
}

double experiment(double* res) {
    double stime;
    double ftime;

    double a1 = 0.0;
    double b1 = 16.0;
    
    double a2 = 0.0;
    double b2 = 16.0;
    
    double h = 0.0001;
    
    stime = clock();
    integral(a1, b1, a2, b2, h, h, res);
    ftime = clock();

    return (stime - ftime) / CLOCKS_PER_SEC;
}

void lab2_start() {
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
