#include <iostream>
#include <iomanip>
#include <vector>
#include <complex>
#include <cmath>
#include <mpi.h>

const int N = 1024; // Число отсчетов
const double PI = 3.14159265358979323846;
const double T = 1.0; // Период

// Функция для генерации тестового сигнала
std::vector<double> generate_signal(double frequency, double amplitude, double duration, int sample_rate)
{
    int samples = static_cast<int>(duration * sample_rate);
    std::vector<double> signal(samples);
    for (int i = 0; i < samples; ++i)
    {
        double t = static_cast<double>(i) / sample_rate;
        signal[i] = amplitude * sin(2 * PI * frequency * t);
    }
    return signal;
}

// Функция для вычисления БПФ
void fft(std::vector<std::complex<double>> &x)
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

    // Рекурсивное вычисление
    fft(even);
    fft(odd);

    // Комбинирование результатов
    for (int i = 0; i < n / 2; ++i)
    {
        double theta = -2.0 * PI * i / n;
        std::complex<double> t = std::polar(1.0, theta) * odd[i];
        x[i] = even[i] + t;
        x[i + n / 2] = even[i] - t;
    }
}

// Нормализация амплитудного спектра
void normalize_spectrum(std::vector<double> &spectrum)
{
    double max_val = *std::max_element(spectrum.begin(), spectrum.end());
    for (auto &val : spectrum)
    {
        val /= max_val;
    }
}

// Вычисление коэффициентов ряда Фурье
std::vector<double> compute_fourier_coefficients(const std::vector<std::complex<double>> &fft_result)
{
    std::vector<double> coefficients(N);
    for (int k = 0; k < N; ++k)
    {
        coefficients[k] = std::abs(fft_result[k]) / N; // Нормируем
    }
    return coefficients;
}

// Вычисление функции с использованием коэффициентов
double calculate_function(double t, const std::vector<double> &coefficients)
{
    double sum = 0.0;
    for (int k = 1; k < coefficients.size(); ++k)
    {
        sum += (coefficients[k] * sin(k * (2 * PI / T) * t));
    }
    return (PI / 2.0) - (PI * t / T) + sum;
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Генерация тестового сигнала
    double frequency = 10.0; // Частота 10 Hz
    double amplitude = 1.0;  // Амплитуда 1
    double duration = 1.0;   // Длительность 1 секунда
    int sample_rate = N;     // Частота дискретизации равна числу отсчетов
    std::vector<double> signal;

    if (rank == 0)
    {
        signal = generate_signal(frequency, amplitude, duration, sample_rate);
    }

    // Распределение сигнала между процессами
    std::vector<double> local_signal(N / size);
    MPI_Scatter(signal.data(), N / size, MPI_DOUBLE, local_signal.data(), N / size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Применение БПФ к локальному сигналу
    std::vector<std::complex<double>> fft_input(local_signal.begin(), local_signal.end());
    fft(fft_input);

    // Собираем результаты БПФ
    std::vector<std::complex<double>> fft_result(N);
    MPI_Gather(fft_input.data(), N / size, MPI_DOUBLE_COMPLEX, fft_result.data(), N / size, MPI_DOUBLE_COMPLEX, 0, MPI_COMM_WORLD);

    // Нормализация спектра
    if (rank == 0)
    {
        std::vector<double> amplitude_spectrum(N);
        for (int i = 0; i < N; ++i)
        {
            amplitude_spectrum[i] = std::abs(fft_result[i]);
        }
        normalize_spectrum(amplitude_spectrum);

        // Вычисление коэффициентов ряда Фурье
        std::vector<double> coefficients = compute_fourier_coefficients(fft_result);

        // Сравнение с аналитическим значением
        std::cout << std::fixed << std::setprecision(6);
        std::cout << "t\t\tCalculated\tAnalytic\n";
        for (double t = 0.1; t < T; t += 0.1)
        {
            double calculated_value = calculate_function(t, coefficients);
            double analytic_value = (PI / 2.0) - (PI * t / T);
            std::cout << t << "\t" << calculated_value << "\t" << analytic_value << "\n";
        }
    }

    MPI_Finalize();
    return 0;
}
