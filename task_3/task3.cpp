#include <omp.h>
#include <iostream>
#include <cmath>
#include <chrono>
#include <random>

#define M_PI 3.14159265358979323846

double f(double x, double y)
{
    return exp(sin(M_PI * x) * cos(M_PI * y)) + 1;
}

double rectangular_integral(double x_min, double x_max, double y_min, double y_max, int nx, int ny)
{
    double dx = (x_max - x_min) / nx;
    double dy = (y_max - y_min) / ny;
    double sum = 0.0;

#pragma omp parallel for reduction(+ : sum)
    for (int i = 0; i < nx; ++i)
    {
        for (int j = 0; j < ny; ++j)
        {
            double x = x_min + (i + 0.5) * dx;
            double y = y_min + (j + 0.5) * dy;
            sum += f(x, y) * dx * dy;
        }
    }

    return sum;
}

double trapezoidal_integral(double x_min, double x_max, double y_min, double y_max, int nx, int ny)
{
    double dx = (x_max - x_min) / nx;
    double dy = (y_max - y_min) / ny;
    double sum = 0.0;

#pragma omp parallel for reduction(+ : sum)
    for (int i = 0; i <= nx; ++i)
    {
        for (int j = 0; j <= ny; ++j)
        {
            double x = x_min + i * dx;
            double y = y_min + j * dy;
            double weight = (i == 0 || i == nx) && (j == 0 || j == ny) ? 1 : (i == 0 || i == nx || j == 0 || j == ny) ? 0.5
                                                                                                                      : 1;
            sum += weight * f(x, y);
        }
    }

    return sum * dx * dy / 4; // Dividing by 4 due to double counting
}

double monte_carlo_integral(double x_min, double x_max, double y_min, double y_max, int samples)
{
    double sum = 0.0;

#pragma omp parallel
    {
        // Создаем генератор случайных чисел для каждого потока
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis_x(x_min, x_max);
        std::uniform_real_distribution<> dis_y(y_min, y_max);

        double local_sum = 0.0;

#pragma omp for
        for (int i = 0; i < samples; ++i)
        {
            double x = dis_x(gen); // Случайное число в диапазоне [x_min, x_max]
            double y = dis_y(gen); // Случайное число в диапазоне [y_min, y_max]
            local_sum += f(x, y);
        }

#pragma omp atomic
        sum += local_sum;
    }

    return sum * (x_max - x_min) * (y_max - y_min) / samples;
}

int main()
{
    double result = rectangular_integral(0, 16, 0, 16, 1000, 1000);
    std::cout << "Rectangular Integral: " << result << std::endl;

    double result = trapezoidal_integral(0, 16, 0, 16, 1000, 1000);
    std::cout << "Trapezoidal Integral: " << result << std::endl;

    double result = monte_carlo_integral(0, 16, 0, 16, 1000000);
    std::cout << "Monte Carlo Integral: " << result << std::endl;

    return 0;
}