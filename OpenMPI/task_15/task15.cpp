#include "includes.h"

const int MAX_ITER = 10000;
const double TOL = 1e-6;

// Последовательный метод Гаусса
void gauss_sequential(std::vector<std::vector<double>> &A, std::vector<double> &b, std::vector<double> &x)
{
    int n = A.size();
    // Прямой ход
    for (int k = 0; k < n; ++k)
    {
        for (int i = k + 1; i < n; ++i)
        {
            double factor = A[i][k] / A[k][k];
            for (int j = k; j < n; ++j)
            {
                A[i][j] -= factor * A[k][j];
            }
            b[i] -= factor * b[k];
        }
    }
    // Обратный ход
    for (int i = n - 1; i >= 0; --i)
    {
        x[i] = b[i];
        for (int j = i + 1; j < n; ++j)
        {
            x[i] -= A[i][j] * x[j];
        }
        x[i] /= A[i][i];
    }
}

// MPI метод Гаусса (распараллеливание по строкам)
void gauss_mpi(std::vector<std::vector<double>> &A, std::vector<double> &b, std::vector<double> &x, int rank, int size)
{
    int n = A.size();
    for (int k = 0; k < n; ++k)
    {
        if (rank == k % size)
        {
            // Широковещательная отправка k-ой строки
            MPI_Bcast(&A[k][0], n, MPI_DOUBLE, k % size, MPI_COMM_WORLD);
            MPI_Bcast(&b[k], 1, MPI_DOUBLE, k % size, MPI_COMM_WORLD);
        }
        else
        {
            MPI_Bcast(&A[k][0], n, MPI_DOUBLE, k % size, MPI_COMM_WORLD);
            MPI_Bcast(&b[k], 1, MPI_DOUBLE, k % size, MPI_COMM_WORLD);
        }
        for (int i = k + 1; i < n; ++i)
        {
            if (rank == i % size)
            {
                double factor = A[i][k] / A[k][k];
                for (int j = k; j < n; ++j)
                {
                    A[i][j] -= factor * A[k][j];
                }
                b[i] -= factor * b[k];
            }
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    for (int i = n - 1; i >= 0; --i)
    {
        if (rank == i % size)
        {
            x[i] = b[i];
            for (int j = i + 1; j < n; ++j)
            {
                x[i] -= A[i][j] * x[j];
            }
            x[i] /= A[i][i];
        }
        MPI_Bcast(&x[i], 1, MPI_DOUBLE, i % size, MPI_COMM_WORLD);
    }
}

// Метод Якоби
void jacobi_sequential(const std::vector<std::vector<double>> &A, const std::vector<double> &b, std::vector<double> &x)
{
    int n = A.size();
    std::vector<double> x_old(n);
    for (int iter = 0; iter < MAX_ITER; ++iter)
    {
        x_old = x;
        for (int i = 0; i < n; ++i)
        {
            double sigma = 0.0;
            for (int j = 0; j < n; ++j)
            {
                if (j != i)
                {
                    sigma += A[i][j] * x_old[j];
                }
            }
            x[i] = (b[i] - sigma) / A[i][i];
        }
        // Проверка на сходимость
        double error = 0.0;
        for (int i = 0; i < n; ++i)
        {
            error += std::pow(x[i] - x_old[i], 2);
        }
        if (std::sqrt(error) < TOL)
        {
            break;
        }
    }
}

// Функция для параллельного метода Якоби с использованием MPI
void jacobi_mpi(const std::vector<std::vector<double>> &A, const std::vector<double> &b, std::vector<double> &x, int rank, int size)
{
    int n = A.size();
    std::vector<double> x_old(n, 0.0);
    std::vector<double> x_local(n / size + (rank < n % size ? 1 : 0), 0.0); // Часть вектора x для каждого процесса

    int rows_per_proc = n / size;                                      // Количество строк на процесс
    int extra_rows = n % size;                                         // Лишние строки при делении нацело
    int start_row = rank * rows_per_proc + std::min(rank, extra_rows); // Начальная строка для текущего процесса
    int local_rows = rows_per_proc + (rank < extra_rows ? 1 : 0);      // Строки, обрабатываемые текущим процессом

    for (int iter = 0; iter < MAX_ITER; ++iter)
    {
        x_old = x; // Сохраняем старые значения для всех процессов
        MPI_Allgather(MPI_IN_PLACE, local_rows, MPI_DOUBLE, &x[0], local_rows, MPI_DOUBLE, MPI_COMM_WORLD);

        for (int i = 0; i < local_rows; ++i)
        {
            int global_i = start_row + i;
            double sigma = 0.0;
            for (int j = 0; j < n; ++j)
            {
                if (j != global_i)
                {
                    sigma += A[global_i][j] * x_old[j];
                }
            }
            x_local[i] = (b[global_i] - sigma) / A[global_i][i];
        }

        // Синхронизация вычисленных значений x с остальными процессами
        MPI_Allgather(&x_local[0], local_rows, MPI_DOUBLE, &x[0], local_rows, MPI_DOUBLE, MPI_COMM_WORLD);

        // Проверка на сходимость (рассчитывается только на процессах)
        double local_error = 0.0;
        for (int i = 0; i < local_rows; ++i)
        {
            local_error += std::pow(x_local[i] - x_old[start_row + i], 2);
        }

        // Глобальная проверка ошибки на всех процессах
        double global_error = 0.0;
        MPI_Allreduce(&local_error, &global_error, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        if (std::sqrt(global_error) < TOL)
        {
            break;
        }
    }
}

// Последовательная версия метода Зейделя
void seidel_sequential(const std::vector<std::vector<double>> &A, const std::vector<double> &b, std::vector<double> &x)
{
    int n = A.size();
    std::vector<double> x_old(n);

    for (int iter = 0; iter < MAX_ITER; ++iter)
    {
        x_old = x;

        for (int i = 0; i < n; ++i)
        {
            double sigma = 0.0;
            for (int j = 0; j < n; ++j)
            {
                if (j != i)
                {
                    sigma += A[i][j] * x[j]; // Используем как старые, так и новые значения x
                }
            }
            x[i] = (b[i] - sigma) / A[i][i];
        }

        // Проверка на сходимость
        double error = 0.0;
        for (int i = 0; i < n; ++i)
        {
            error += std::pow(x[i] - x_old[i], 2);
        }

        if (std::sqrt(error) < TOL)
        {
            break;
        }
    }
}

// MPI метод Зейделя
void seidel_mpi(const std::vector<std::vector<double>> &A, const std::vector<double> &b, std::vector<double> &x, int rank, int size)
{
    int n = A.size();
    std::vector<double> x_old(n, 0.0);
    std::vector<double> x_local(n / size + (rank < n % size ? 1 : 0), 0.0); // Локальный кусок x для каждого процесса

    int rows_per_proc = n / size;                                      // Количество строк на процесс
    int extra_rows = n % size;                                         // Лишние строки
    int start_row = rank * rows_per_proc + std::min(rank, extra_rows); // Начальная строка для процесса
    int local_rows = rows_per_proc + (rank < extra_rows ? 1 : 0);      // Количество строк для текущего процесса

    for (int iter = 0; iter < MAX_ITER; ++iter)
    {
        x_old = x; // Сохраняем старые значения для всех процессов

        for (int i = 0; i < local_rows; ++i)
        {
            int global_i = start_row + i;
            double sigma = 0.0;

            // Расчет суммы по строкам (используем частично новые и старые значения)
            for (int j = 0; j < n; ++j)
            {
                if (j != global_i)
                {
                    sigma += A[global_i][j] * x[j]; // Использование глобального вектора x
                }
            }
            x_local[i] = (b[global_i] - sigma) / A[global_i][i];

            // Синхронизация обновленного значения с остальными процессами
            x[global_i] = x_local[i]; // Немедленное обновление глобального вектора x
            MPI_Allgather(MPI_IN_PLACE, local_rows, MPI_DOUBLE, &x[0], local_rows, MPI_DOUBLE, MPI_COMM_WORLD);
        }

        // Проверка на сходимость (рассчитывается только на процессах)
        double local_error = 0.0;
        for (int i = 0; i < local_rows; ++i)
        {
            local_error += std::pow(x_local[i] - x_old[start_row + i], 2);
        }

        // Глобальная проверка ошибки
        double global_error = 0.0;
        MPI_Allreduce(&local_error, &global_error, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        if (std::sqrt(global_error) < TOL)
        {
            break;
        }
    }
}

// Последовательная версия метода SOR
void sor_sequential(const std::vector<std::vector<double>> &A, const std::vector<double> &b, std::vector<double> &x, double omega)
{
    int n = A.size();
    std::vector<double> x_old(n);

    for (int iter = 0; iter < MAX_ITER; ++iter)
    {
        x_old = x;

        for (int i = 0; i < n; ++i)
        {
            double sigma = 0.0;
            for (int j = 0; j < n; ++j)
            {
                if (j != i)
                {
                    sigma += A[i][j] * x[j]; // Используем старые значения
                }
            }
            x[i] = (1.0 - omega) * x_old[i] + (omega * (b[i] - sigma) / A[i][i]);
        }

        // Проверка на сходимость
        double error = 0.0;
        for (int i = 0; i < n; ++i)
        {
            error += std::pow(x[i] - x_old[i], 2);
        }

        if (std::sqrt(error) < TOL)
        {
            break;
        }
    }
}

// MPI метод верхней релаксации (SOR)
void sor_mpi(const std::vector<std::vector<double>> &A, const std::vector<double> &b, std::vector<double> &x, double omega, int rank, int size)
{
    int n = A.size();
    std::vector<double> x_old(n, 0.0);
    std::vector<double> x_local(n / size + (rank < n % size ? 1 : 0), 0.0);

    int rows_per_proc = n / size;                                      // Количество строк на процесс
    int extra_rows = n % size;                                         // Лишние строки
    int start_row = rank * rows_per_proc + std::min(rank, extra_rows); // Начальная строка
    int local_rows = rows_per_proc + (rank < extra_rows ? 1 : 0);      // Количество строк для текущего процесса

    for (int iter = 0; iter < MAX_ITER; ++iter)
    {
        x_old = x; // Сохраняем старые значения для всех процессов

        for (int i = 0; i < local_rows; ++i)
        {
            int global_i = start_row + i;
            double sigma = 0.0;

            // Расчет суммы по строкам (используем старые значения)
            for (int j = 0; j < n; ++j)
            {
                if (j != global_i)
                {
                    sigma += A[global_i][j] * x[j]; // Используем глобальный вектор x
                }
            }
            x_local[i] = (1.0 - omega) * x_old[global_i] + (omega * (b[global_i] - sigma) / A[global_i][global_i]);

            // Синхронизация обновленного значения с остальными процессами
            x[global_i] = x_local[i]; // Обновляем глобальный вектор x
            MPI_Allgather(MPI_IN_PLACE, local_rows, MPI_DOUBLE, &x[0], local_rows, MPI_DOUBLE, MPI_COMM_WORLD);
        }

        // Проверка на сходимость (рассчитывается только на процессах)
        double local_error = 0.0;
        for (int i = 0; i < local_rows; ++i)
        {
            local_error += std::pow(x_local[i] - x_old[start_row + i], 2);
        }

        // Глобальная проверка ошибки
        double global_error = 0.0;
        MPI_Allreduce(&local_error, &global_error, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        if (std::sqrt(global_error) < TOL)
        {
            break;
        }
    }
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Пример настройки задачи
    int n = 100;         // Размерность системы
    double omega = 1.25; // Параметр релаксации для SOR
    std::vector<std::vector<double>> A(n, std::vector<double>(n));
    std::vector<double> b(n), x(n, 0.0); // Вектор решений и правых частей

    // Инициализация матрицы A и вектора b
    // ...

    // Таблица 1
    if (rank == 0)
    {
        std::cout << std::fixed << std::setprecision(6);
        std::cout << "Таблица 1\n";
        std::cout << "---------------------------------------------------------------------------------------------------------\n";
        std::cout << "| Номер теста | Порядок системы | Последовательный алгоритм | Параллельный алгоритм | Время | Ускорение |\n";
        std::cout << "---------------------------------------------------------------------------------------------------------\n";
    }

    // Последовательные методы
    auto start = MPI_Wtime();
    gauss_sequential(A, b, x);
    double gauss_time_seq = MPI_Wtime() - start;

    start = MPI_Wtime();
    jacobi_sequential(A, b, x);
    double jacobi_time_seq = MPI_Wtime() - start;

    start = MPI_Wtime();
    seidel_sequential(A, b, x);
    double seidel_time_seq = MPI_Wtime() - start;

    start = MPI_Wtime();
    sor_sequential(A, b, x, omega);
    double sor_time_seq = MPI_Wtime() - start;

    // Параллельные методы
    start = MPI_Wtime();
    gauss_mpi(A, b, x, rank, size);
    double gauss_time_mpi = MPI_Wtime() - start;

    start = MPI_Wtime();
    jacobi_mpi(A, b, x, rank, size);
    double jacobi_time_mpi = MPI_Wtime() - start;

    start = MPI_Wtime();
    seidel_mpi(A, b, x, rank, size);
    double seidel_time_mpi = MPI_Wtime() - start;

    start = MPI_Wtime();
    sor_mpi(A, b, x, omega, rank, size);
    double sor_time_mpi = MPI_Wtime() - start;

    // Собираем результаты на процессе 0
    if (rank == 0)
    {
        std::cout << "1           | " << n << "             | " << gauss_time_seq << "           | " << gauss_time_mpi << "           | "
                  << gauss_time_seq / gauss_time_mpi << "\n";
        std::cout << "2           | " << n << "             | " << jacobi_time_seq << "           | " << jacobi_time_mpi << "           | "
                  << jacobi_time_seq / jacobi_time_mpi << "\n";
        std::cout << "3           | " << n << "             | " << seidel_time_seq << "           | " << seidel_time_mpi << "           | "
                  << seidel_time_seq / seidel_time_mpi << "\n";
        std::cout << "4           | " << n << "             | " << sor_time_seq << "           | " << sor_time_mpi << "           | "
                  << sor_time_seq / sor_time_mpi << "\n";
    }

    // Таблица 2 (ускорение)
    if (rank == 0)
    {
        std::cout << "\nТаблица 2\n";
        std::cout << "----------------------------------------------------------------------------------------------------------\n";
        std::cout << "| Номер теста | Порядок системы | Ускорение Гаусса | Ускорение Якоби | Ускорение Зейделя | Ускорение SOR |\n";
        std::cout << "----------------------------------------------------------------------------------------------------------\n";
        std::cout << "1           | " << n << "             | " << gauss_time_seq / gauss_time_mpi << "           | "
                  << jacobi_time_seq / jacobi_time_mpi << "           | " << seidel_time_seq / seidel_time_mpi << "           | "
                  << sor_time_seq / sor_time_mpi << "\n";
    }

    MPI_Finalize();
    return 0;
}
