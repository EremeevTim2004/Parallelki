#include "includes.h"

// Константы для интегрирования
const double a1 = 0.0, b1 = 16.0; // Пределы интегрирования первого интеграла
const double a2 = 0.0, b2 = 16.0; // Пределы интегрирования второго интеграла
const int N = 1000000;            // Количество точек или разбиений

// Тестовая функция для интегрирования
double func(double x, double y)
{
    return (exp(sin(M_PI * x) * cos(M_PI * y)) + 1) / ((b1 - a1) * (b2 - a2));
}

// Метод ячеек
double cell_method(double a1, double b1, double a2, double b2, int n, int rank, int size)
{
    double h1 = (b1 - a1) / n; // Шаг по x
    double h2 = (b2 - a2) / n; // Шаг по y
    double local_sum = 0.0;

    int local_n = n / size;
    double local_a1 = a1 + rank * local_n * h1;
    double local_b1 = local_a1 + local_n * h1;

    for (int i = 0; i < local_n; ++i)
    {
        double x = local_a1 + i * h1 + h1 / 2.0;
        for (int j = 0; j < n; ++j)
        {
            double y = a2 + j * h2 + h2 / 2.0;
            local_sum += func(x, y);
        }
    }

    return local_sum * h1 * h2;
}

// Метод трапеций
double trapezoidal_method(double a1, double b1, double a2, double b2, int n, int rank, int size)
{
    double h1 = (b1 - a1) / n; // Шаг по x
    double h2 = (b2 - a2) / n; // Шаг по y
    double local_sum = 0.0;

    int local_n = n / size;
    double local_a1 = a1 + rank * local_n * h1;
    double local_b1 = local_a1 + local_n * h1;

    for (int i = 0; i <= local_n; ++i)
    {
        double x = local_a1 + i * h1;
        for (int j = 0; j <= n; ++j)
        {
            double y = a2 + j * h2;
            double weight = ((i == 0 || i == local_n) ? 0.5 : 1.0) * ((j == 0 || j == n) ? 0.5 : 1.0);
            local_sum += weight * func(x, y);
        }
    }

    return local_sum * h1 * h2;
}

// Метод статистических испытаний (Монте-Карло)
double monte_carlo_method(double a1, double b1, double a2, double b2, int n, int rank, int size)
{
    double local_sum = 0.0;

    // Установить зерно для случайного числа
    srand(time(NULL) + rank);

    int local_n = n / size;

    for (int i = 0; i < local_n; ++i)
    {
        double x = a1 + (b1 - a1) * (rand() / (double)RAND_MAX);
        double y = a2 + (b2 - a2) * (rand() / (double)RAND_MAX);
        local_sum += func(x, y);
    }

    return local_sum * (b1 - a1) * (b2 - a2) / n;
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int n = N; // Количество точек или разбиений

    double local_result = 0.0, global_result = 0.0;
    double start_time, end_time;

    // 1. Метод ячеек
    if (rank == 0)
        std::cout << "Метод ячеек:\n";
    start_time = MPI_Wtime();
    local_result = cell_method(a1, b1, a2, b2, n, rank, size);
    MPI_Reduce(&local_result, &global_result, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    end_time = MPI_Wtime();
    if (rank == 0)
    {
        std::cout << "Численное значение: " << global_result << std::endl;
        std::cout << "Время выполнения: " << end_time - start_time << " секунд.\n";
    }

    // 2. Метод трапеций
    if (rank == 0)
        std::cout << "\nМетод трапеций:\n";
    start_time = MPI_Wtime();
    local_result = trapezoidal_method(a1, b1, a2, b2, n, rank, size);
    MPI_Reduce(&local_result, &global_result, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    end_time = MPI_Wtime();
    if (rank == 0)
    {
        std::cout << "Численное значение: " << global_result << std::endl;
        std::cout << "Время выполнения: " << end_time - start_time << " секунд.\n";
    }

    // 3. Метод Монте-Карло
    if (rank == 0)
        std::cout << "\nМетод Монте-Карло:\n";
    start_time = MPI_Wtime();
    local_result = monte_carlo_method(a1, b1, a2, b2, n, rank, size);
    MPI_Reduce(&local_result, &global_result, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    end_time = MPI_Wtime();
    if (rank == 0)
    {
        std::cout << "Численное значение: " << global_result << std::endl;
        std::cout << "Время выполнения: " << end_time - start_time << " секунд.\n";
    }

    MPI_Finalize();
    return 0;
}