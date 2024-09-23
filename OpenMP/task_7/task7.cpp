#include "includes.h"

// Последовательное БПФ
void fft_recursive(vector<complex<double>> &x)
{
    int n = x.size();
    if (n <= 1)
        return;

    vector<complex<double>> even(n / 2), odd(n / 2);
    for (int i = 0; i < n / 2; ++i)
    {
        even[i] = x[i * 2];
        odd[i] = x[i * 2 + 1];
    }

    fft_recursive(even);
    fft_recursive(odd);

    for (int k = 0; k < n / 2; ++k)
    {
        complex<double> t = polar(1.0, -2 * M_PI * k / n) * odd[k];
        x[k] = even[k] + t;
        x[k + n / 2] = even[k] - t;
    }
}

// Параллельное БПФ
void fft_parallel(vector<complex<double>> &x)
{
    int n = x.size();
    if (n <= 1)
        return;

    vector<complex<double>> even(n / 2), odd(n / 2);
#pragma omp parallel
    {
#pragma omp for
        for (int i = 0; i < n / 2; ++i)
        {
            even[i] = x[i * 2];
            odd[i] = x[i * 2 + 1];
        }
    }

    fft_parallel(even);
    fft_parallel(odd);

#pragma omp parallel for
    for (int k = 0; k < n / 2; ++k)
    {
        complex<double> t = polar(1.0, -2 * M_PI * k / n) * odd[k];
        x[k] = even[k] + t;
        x[k + n / 2] = even[k] - t;
    }
}

int main()
{
    vector<int> sizes = {32768, 65536, 131072, 262144, 524288};
    cout << "Тест\tРазмер входного сигнала\tМин. время (сек)\tМин. время паралл.\tУскорение\n";

    for (size_t i = 0; i < sizes.size(); ++i)
    {
        int size = sizes[i];
        vector<complex<double>> data(size);

        // Заполнение случайными данными
        for (int j = 0; j < size; ++j)
        {
            data[j] = complex<double>(rand() % 100, rand() % 100);
        }

        // Последовательное выполнение
        auto start = chrono::high_resolution_clock::now();
        fft_recursive(data);
        auto end = chrono::high_resolution_clock::now();
        double seq_time = chrono::duration<double>(end - start).count();

        // Параллельное выполнение
        for (int j = 0; j < size; ++j)
        {
            data[j] = complex<double>(rand() % 100, rand() % 100); // Снова заполняем данными
        }

        start = chrono::high_resolution_clock::now();
        fft_parallel(data);
        end = chrono::high_resolution_clock::now();
        double par_time = chrono::duration<double>(end - start).count();

        // Вычисление ускорения
        double speedup = seq_time / par_time;

        // Вывод результатов
        cout << i + 1 << "\t" << size << "\t\t\t"
             << seq_time << "\t\t" << par_time << "\t\t" << speedup << "\n";
    }

    return 0;
}