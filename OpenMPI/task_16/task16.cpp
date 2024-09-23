#include "includes.h"

const int MAX_SIZE = 1048576; // Максимальный размер сигнала
const int TESTS = 5;          // Количество тестов
const double PI = 3.14159265358979323846;

void fft_sequential(std::vector<std::complex<double>> &x)
{
    int n = x.size();
    if (n <= 1)
        return;

    // Разделяем на четные и нечетные элементы
    std::vector<std::complex<double>> even(n / 2), odd(n / 2);
    for (int i = 0; i < n / 2; ++i)
    {
        even[i] = x[i * 2];
        odd[i] = x[i * 2 + 1];
    }

    // Рекурсивное преобразование
    fft_sequential(even);
    fft_sequential(odd);

    // Комбинирование результатов
    for (int i = 0; i < n / 2; ++i)
    {
        double theta = -2.0 * PI * i / n;
        std::complex<double> t = std::polar(1.0, theta) * odd[i];
        x[i] = even[i] + t;
        x[i + n / 2] = even[i] - t;
    }
}

void fft_parallel(std::vector<std::complex<double>> &x, int rank, int size)
{
    int n = x.size();
    int local_n = n / size;
    std::vector<std::complex<double>> local_x(local_n);

    // Распределение данных
    MPI_Scatter(x.data(), local_n * sizeof(std::complex<double>), MPI_BYTE,
                local_x.data(), local_n * sizeof(std::complex<double>), MPI_BYTE, 0, MPI_COMM_WORLD);

    // Локальное вычисление БПФ
    fft_sequential(local_x);

    // Собираем результаты
    MPI_Gather(local_x.data(), local_n * sizeof(std::complex<double>), MPI_BYTE,
               x.data(), local_n * sizeof(std::complex<double>), MPI_BYTE, 0, MPI_COMM_WORLD);
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    std::vector<int> sizes = {1024, 4096, 16384, 65536, 1048576}; // Размеры входного сигнала
    std::vector<double> seq_times(TESTS), par_times(TESTS);

    if (rank == 0)
    {
        std::cout << std::fixed << std::setprecision(6);
        std::cout << "Таблица результатов БПФ\n";
        std::cout << "---------------------------------------------------------------------------------------------------\n";
        std::cout << "| Номер теста | Размер входного сигнала | Мин. время (сек) | Параллельное время (сек) | Ускорение |\n";
        std::cout << "---------------------------------------------------------------------------------------------------\n";
    }

    for (int test = 0; test < TESTS; ++test)
    {
        int n = sizes[test];
        std::vector<std::complex<double>> signal(n);

        // Инициализация входного сигнала
        for (int i = 0; i < n; ++i)
        {
            signal[i] = std::complex<double>(static_cast<double>(i), 0);
        }

        // Последовательное выполнение
        double start_time = MPI_Wtime();
        fft_sequential(signal);
        double end_time = MPI_Wtime();
        seq_times[test] = end_time - start_time;

        // Параллельное выполнение
        start_time = MPI_Wtime();
        fft_parallel(signal, rank, size);
        end_time = MPI_Wtime();
        par_times[test] = end_time - start_time;

        // Вывод результатов на процессе 0
        if (rank == 0)
        {
            double speedup = seq_times[test] / par_times[test];
            std::cout << test + 1 << "           | " << n << "                  | "
                      << seq_times[test] << "           | " << par_times[test] << "             | "
                      << speedup << "\n";
        }
    }

    MPI_Finalize();
    return 0;
}
