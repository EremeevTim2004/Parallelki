#include "includes.h"

// Константы для интегрирования
const double a1 = 0.0, b1 = 16.0;
const double a2 = 0.0, b2 = 16.0;
const int N = 1000000; // Количество отрезков разбиения для каждого измерения

// Тестовая функция для интегрирования
double func(double x, double y)
{
    return (exp(sin(M_PI * x) * cos(M_PI * y)) + 1) / ((b1 - a1) * (b2 - a2));
}

// Метод прямоугольников
double rectangle_method(double a1, double b1, double a2, double b2, int n, int rank, int size)
{
    double h1 = (b1 - a1) / n; // Шаг по x
    double h2 = (b2 - a2) / n; // Шаг по y
    double local_sum = 0.0;

    // Количество отрезков для каждого процесса
    int local_n = n / size;

    // Локальные границы для каждого процесса
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

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int n = N; // Число разбиений
    double local_result = 0.0, global_result = 0.0;

    // Измерение времени
    auto start = std::chrono::high_resolution_clock::now();

    // Вызов метода прямоугольников
    local_result = rectangle_method(a1, b1, a2, b2, n, rank, size);

    // Сбор результатов на главном процессе
    MPI_Reduce(&local_result, &global_result, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    // Измерение времени окончания
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;

    // Вывод результата и времени на главном процессе
    if (rank == 0)
    {
        std::cout << "Численное значение интеграла: " << global_result << std::endl;
        std::cout << "Время выполнения: " << duration.count() << " секунд." << std::endl;
    }

    MPI_Finalize();
    return 0;
}