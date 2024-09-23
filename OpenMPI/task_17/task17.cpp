#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>
#include <mpi.h>

// Функция для бит-реверсирования
unsigned int bit_reverse(unsigned int num, int bits)
{
    unsigned int reversed = 0;
    for (int i = 0; i < bits; ++i)
    {
        if (num & (1 << i))
        {
            reversed |= (1 << ((bits - 1) - i));
        }
    }
    return reversed;
}

void bit_reverse_parallel(std::vector<unsigned int> &data, int rank, int size)
{
    int n = data.size();
    int local_n = n / size; // Количество элементов на процесс

    // Локальный массив для хранения результатов
    std::vector<unsigned int> local_data(local_n);

    // Распределение данных между процессами
    MPI_Scatter(data.data(), local_n, MPI_UNSIGNED, local_data.data(), local_n, MPI_UNSIGNED, 0, MPI_COMM_WORLD);

    // Бит-реверсирование локальных данных
    for (int i = 0; i < local_n; ++i)
    {
        local_data[i] = bit_reverse(local_data[i], (int)log2(n));
    }

    // Сбор результатов
    MPI_Gather(local_data.data(), local_n, MPI_UNSIGNED, data.data(), local_n, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    const int N = 1024; // Размер данных (должен быть степенью двойки)
    std::vector<unsigned int> data(N);

    // Инициализация данных (например, от 0 до N-1)
    if (rank == 0)
    {
        for (int i = 0; i < N; ++i)
        {
            data[i] = i;
        }
    }

    // Замер времени выполнения
    double start_time = MPI_Wtime();

    // Вызов функции бит-реверсирования
    bit_reverse_parallel(data, rank, size);

    double end_time = MPI_Wtime();

    // Вывод результатов на процессе 0
    if (rank == 0)
    {
        std::cout << std::fixed << std::setprecision(6);
        std::cout << "Время выполнения параллельного бит-реверсирования: " << (end_time - start_time) << " секунд.\n";
        std::cout << "Результат (первые 10 элементов): ";
        for (int i = 0; i < 10; ++i)
        {
            std::cout << data[i] << " ";
        }
        std::cout << "\n";
    }

    MPI_Finalize();
    return 0;
}
