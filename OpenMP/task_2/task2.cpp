#include <omp.h>
#include <iostream>
#include <cmath>
#include <chrono>

#define M_PI 3.14159265358979323846

double f(double x, double y, double a1, double b1, double a2, double b2)
{
    return (exp(sin(M_PI * x) * cos(M_PI * y)) + 1) / ((b1 - a1) * (b2 - a2));
}

double rectangle_method_2d(double a1, double b1, double a2, double b2, int n)
{
    double hx = (b1 - a1) / n;
    double hy = (b2 - a2) / n;
    double sum = 0.0;

    for (int i = 0; i < n; ++i)
    {
        double x = a1 + i * hx;
        for (int j = 0; j < n; ++j)
        {
            double y = a2 + j * hy;
            sum += f(x, y, a1, b1, a2, b2);
        }
    }

    return sum * hx * hy;
}

double parallel_rectangle_method_2d_x(double a1, double b1, double a2, double b2, int n)
{
    double hx = (b1 - a1) / n;
    double hy = (b2 - a2) / n;
    double sum = 0.0;

#pragma omp parallel for reduction(+ : sum)
    for (int i = 0; i < n; ++i)
    {
        double x = a1 + i * hx;
        for (int j = 0; j < n; ++j)
        {
            double y = a2 + j * hy;
            sum += f(x, y, a1, b1, a2, b2);
        }
    }

    return sum * hx * hy;
}

double parallel_rectangle_method_2d_xy(double a1, double b1, double a2, double b2, int n)
{
    double hx = (b1 - a1) / n;
    double hy = (b2 - a2) / n;
    double sum = 0.0;

#pragma omp parallel for collapse(2) reduction(+ : sum)
    for (int i = 0; i < n; ++i)
    {
        for (int j = 0; j < n; ++j)
        {
            double x = a1 + i * hx;
            double y = a2 + j * hy;
            sum += f(x, y, a1, b1, a2, b2);
        }
    }

    return sum * hx * hy;
}

int main()
{
    double a1 = 0.0, b1 = 16.0;
    double a2 = 0.0, b2 = 16.0;
    int n = 1000; // Количество разбиений

    auto start = std::chrono::high_resolution_clock::now();
    double result = rectangle_method_2d(a1, b1, a2, b2, n);
    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Sequential Rectangle Method 2D: " << result << std::endl;
    std::cout << "Time: " << elapsed.count() << " seconds" << std::endl;

    auto start = std::chrono::high_resolution_clock::now();
    double result = parallel_rectangle_method_2d_x(a1, b1, a2, b2, n);
    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Parallel Rectangle Method 2D (x-axis): " << result << std::endl;
    std::cout << "Time: " << elapsed.count() << " seconds" << std::endl;

    auto start = std::chrono::high_resolution_clock::now();
    double result = parallel_rectangle_method_2d_xy(a1, b1, a2, b2, n);
    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Parallel Rectangle Method 2D (x and y axes): " << result << std::endl;
    std::cout << "Time: " << elapsed.count() << " seconds" << std::endl;

    return 0;
}
