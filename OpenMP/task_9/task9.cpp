#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <fftw3.h>

#define M_PI 3.14159265358979323846

const int N = 1024;     // Number of samples
const double fs = 1024; // Sampling frequency (Hz)
const double T = 1.0;   // Period (s)

double signal_function(double t)
{
    return M_PI / 2 - (M_PI * t / T);
}

int main()
{
    std::vector<double> signal(N);

    // Создание сигнала
    for (int n = 0; n < N; ++n)
    {
        double t = n / fs;
        signal[n] = signal_function(t);
    }

    // Расчет БПФ
    std::vector<std::complex<double>> fft_result(N);
    fftw_plan plan = fftw_plan_dft_r2c_1d(N, signal.data(), reinterpret_cast<fftw_complex *>(fft_result.data()), FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);

    // Расчет коэффициентов ряда Фурье
    std::vector<double> fourier_coefficients(N / 2);
    for (int k = 0; k < N / 2; ++k)
    {
        fourier_coefficients[k] = std::abs(fft_result[k]) / N;
    }

    // Сравнение с аналитическим разложением
    std::cout << "Сравнение с аналитическим разложением:\n";
    for (int k = 1; k <= 10; ++k)
    {                                      // Первые 10 коэффициентов
        double analytical_value = 1.0 / k; // Аналитическое значение
        std::cout << "k: " << k << ", Коэффициент ряда Фурье: " << fourier_coefficients[k - 1] << ", Аналитическое: " << analytical_value << std::endl;
    }

    return 0;
}
