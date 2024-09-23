#include "includes.h"

const int N = 128;           // Number of samples
const double fs = N;         // Sampling frequency (Hz)
const double f = 10.0;       // Signal frequency (Hz)
const double duration = 1.0; // Signal duration (s)

int main()
{
    std::vector<double> signal(N);

    // Создание сигнала
    for (int n = 0; n < N; ++n)
    {
        double t = n / fs;
        signal[n] = sin(2 * M_PI * f * t);
    }

    // Расчет БПФ
    std::vector<std::complex<double>> fft_result(N);
    fftw_plan plan = fftw_plan_dft_r2c_1d(N, signal.data(), reinterpret_cast<fftw_complex *>(fft_result.data()), FFTW_ESTIMATE);
    fftw_execute(plan);
    fftw_destroy_plan(plan);

    // Нормализация амплитудного спектра
    std::vector<double> amplitude_spectrum(N / 2);
    for (int k = 0; k < N / 2; ++k)
    {
        amplitude_spectrum[k] = std::abs(fft_result[k]) / N;
    }

    // Вывод амплитудного спектра
    for (int k = 0; k < N / 2; ++k)
    {
        std::cout << "Frequency: " << k * (fs / N) << " Hz, Amplitude: " << amplitude_spectrum[k] << std::endl;
    }

    return 0;
}
