#include <iostream>
#include <cmath>
#include <mpi.h>

const double infinity = 1e6; // Замена для бесконечности
const int N = 1000000;       // Количество отрезков разбиения

// Функция для интегрирования
double f(double x)
{
    return log(x); // Используем логарифм
}

// Метод прямоугольников
double rectangle_method(double a, double b, int n, int rank, int size)
{
    double h = (b - a) / n; // Шаг разбиения
    double local_sum = 0.0;

    // Определяем локальные границы для каждого процесса
    int local_n = n / size;
    double local_a = a + rank * local_n * h;
    double local_b = local_a + local_n * h;

    for (int i = 0; i < local_n; ++i)
    {
        double x = local_a + i * h + h / 2.0;
        local_sum += f(x);
    }
    return local_sum * h;
}

// Формула Симпсона
double simpson_method(double a, double b, int n, int rank, int size)
{
    if (n % 2 != 0)
    {
        if (rank == 0)
            std::cerr << "n должно быть четным для метода Симпсона.\n";
        return 0.0;
    }

    double h = (b - a) / n;
    double local_sum = 0.0;

    // Определяем локальные границы для каждого процесса
    int local_n = n / size;
    double local_a = a + rank * local_n * h;
    double local_b = local_a + local_n * h;

    local_sum += f(local_a) + f(local_b); // Добавляем крайние точки

    for (int i = 1; i < local_n; ++i)
    {
        double x = local_a + i * h;
        local_sum += (i % 2 == 0) ? 2 * f(x) : 4 * f(x);
    }

    return local_sum * h / 3.0;
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    double a = 0.0;      // Нижний предел интегрирования
    double b = infinity; // Верхний предел (замена для бесконечности)
    int n = N;

    // Локальные результаты для каждого процесса
    double local_rect = rectangle_method(a, b, n, rank, size);
    double local_simp = simpson_method(a, b, n, rank, size);

    // Глобальные результаты на процессе 0
    double global_rect = 0.0;
    double global_simp = 0.0;

    MPI_Reduce(&local_rect, &global_rect, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&local_simp, &global_simp, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0)
    {
        std::cout << "Метод прямоугольников: " << global_rect << std::endl;
        std::cout << "Метод Симпсона: " << global_simp << std::endl;
    }

    MPI_Finalize();

    return 0;
}
