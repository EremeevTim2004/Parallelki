#include "includes.h"

double f(double x)
{
    return log(x); // Используем логарифм
}

double rectangle_method(double a, double b, int n)
{
    double h = (b - a) / n;
    double sum = 0.0;
    for (int i = 0; i < n; ++i)
    {
        sum += f(a + i * h);
    }
    return sum * h;
}

double parallel_rectangle_method(double a, double b, int n)
{
    double h = (b - a) / n;
    double sum = 0.0;

#pragma omp parallel for reduction(+ : sum)
    for (int i = 0; i < n; ++i)
    {
        sum += f(a + i * h);
    }

    return sum * h;
}

double simpson_method(double a, double b, int n)
{
    double h = (b - a) / n;
    double sum = f(a) + f(b);

    for (int i = 1; i < n; i += 2)
    {
        sum += 4 * f(a + i * h);
    }

    for (int i = 2; i < n - 1; i += 2)
    {
        sum += 2 * f(a + i * h);
    }

    return sum * h / 3.0;
}

double parallel_simpson_method(double a, double b, int n)
{
    double h = (b - a) / n;
    double sum = f(a) + f(b);

#pragma omp parallel for reduction(+ : sum)
    for (int i = 1; i < n; i += 2)
    {
        sum += 4 * f(a + i * h);
    }

#pragma omp parallel for reduction(+ : sum)
    for (int i = 2; i < n - 1; i += 2)
    {
        sum += 2 * f(a + i * h);
    }

    return sum * h / 3.0;
}

int main()
{
    double a = 1e-6; // Избегаем log(0), берем очень маленькое значение вместо 0
    double b = 1.0;
    int n = 1000000; // Количество разбиений

    auto start = std::chrono::high_resolution_clock::now();
    double result = rectangle_method(a, b, n);
    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Sequential Rectangle Method: " << result << std::endl;
    std::cout << "Time: " << elapsed.count() << " seconds" << std::endl;

    auto start = std::chrono::high_resolution_clock::now();
    double result = parallel_rectangle_method(a, b, n);
    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Parallel Rectangle Method: " << result << std::endl;
    std::cout << "Time: " << elapsed.count() << " seconds" << std::endl;

    auto start = std::chrono::high_resolution_clock::now();
    double result = simpson_method(a, b, n);
    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Sequential Simpson's Method: " << result << std::endl;
    std::cout << "Time: " << elapsed.count() << " seconds" << std::endl;

    auto start = std::chrono::high_resolution_clock::now();
    double result = parallel_simpson_method(a, b, n);
    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Parallel Simpson's Method: " << result << std::endl;
    std::cout << "Time: " << elapsed.count() << " seconds" << std::endl;

    return 0;
}
