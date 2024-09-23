#include "includes.h"

const double TOL = 1e-6;   // Заданная точность
const int MAX_ITER = 1000; // Максимальное количество итераций

// Прямой метод Гаусса для приведения к верхнетреугольному виду
void gauss_forward_elimination(std::vector<std::vector<double>> &A, std::vector<double> &b)
{
    int n = A.size();

    for (int k = 0; k < n; ++k)
    {
        // Поиск строки с максимальным элементом для устойчивости метода
        int maxRow = k;
        for (int i = k + 1; i < n; ++i)
        {
            if (std::abs(A[i][k]) > std::abs(A[maxRow][k]))
            {
                maxRow = i;
            }
        }

        // Меняем строки местами
        std::swap(A[k], A[maxRow]);
        std::swap(b[k], b[maxRow]);

// Параллельное обнуление элементов ниже диагонали в k-м столбце
#pragma omp parallel for
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
}

// Обратная подстановка для нахождения решения системы
void gauss_back_substitution(const std::vector<std::vector<double>> &A, const std::vector<double> &b, std::vector<double> &x)
{
    int n = A.size();

    for (int i = n - 1; i >= 0; --i)
    {
        x[i] = b[i];
#pragma omp parallel for reduction(- : x[i])
        for (int j = i + 1; j < n; ++j)
        {
            x[i] -= A[i][j] * x[j];
        }
        x[i] /= A[i][i];
    }
}

void jacobi_parallel(const std::vector<std::vector<double>> &A, const std::vector<double> &b, std::vector<double> &x)
{
    int n = A.size();
    std::vector<double> x_new(n, 0.0); // Вектор для хранения новых значений переменных

    int iter = 0;
    double diff = 0.0;

    do
    {
        diff = 0.0;

#pragma omp parallel for reduction(max : diff)
        for (int i = 0; i < n; ++i)
        {
            double sum = 0.0;
            for (int j = 0; j < n; ++j)
            {
                if (i != j)
                {
                    sum += A[i][j] * x[j];
                }
            }

            // Обновляем значение x_new для i-й переменной
            x_new[i] = (b[i] - sum) / A[i][i];

            // Находим максимальное изменение для оценки сходимости
            diff = std::max(diff, std::abs(x_new[i] - x[i]));
        }

        // Копируем новые значения в вектор x для следующей итерации
        x = x_new;
        ++iter;
    } while (diff > TOL && iter < MAX_ITER);

    if (iter == MAX_ITER)
    {
        std::cout << "Метод Якоби не сошелся за " << MAX_ITER << " итераций.\n";
    }
    else
    {
        std::cout << "Метод Якоби сошелся за " << iter << " итераций.\n";
    }
}

void gauss_seidel_parallel(const std::vector<std::vector<double>> &A, const std::vector<double> &b, std::vector<double> &x)
{
    int n = A.size();
    int iter = 0;
    double diff = 0.0;

    do
    {
        diff = 0.0;

#pragma omp parallel for reduction(max : diff)
        for (int i = 0; i < n; ++i)
        {
            double sum = 0.0;
            for (int j = 0; j < i; ++j)
            {
                sum += A[i][j] * x[j]; // Используем обновленные значения
            }
            for (int j = i + 1; j < n; ++j)
            {
                sum += A[i][j] * x[j]; // Используем старые значения
            }

            double new_xi = (b[i] - sum) / A[i][i];
            diff = std::max(diff, std::abs(new_xi - x[i]));

#pragma omp critical
            x[i] = new_xi; // Обновляем переменную
        }

        ++iter;
    } while (diff > TOL && iter < MAX_ITER);

    if (iter == MAX_ITER)
    {
        std::cout << "Метод Зейделя не сошелся за " << MAX_ITER << " итераций.\n";
    }
    else
    {
        std::cout << "Метод Зейделя сошелся за " << iter << " итераций.\n";
    }
}

void sor_parallel(const std::vector<std::vector<double>> &A, const std::vector<double> &b, std::vector<double> &x, double omega)
{
    int n = A.size();
    int iter = 0;
    double diff = 0.0;

    do
    {
        diff = 0.0;

#pragma omp parallel for reduction(max : diff)
        for (int i = 0; i < n; ++i)
        {
            double sum1 = 0.0;
            for (int j = 0; j < i; ++j)
            {
                sum1 += A[i][j] * x[j]; // Используем уже обновленные значения
            }

            double sum2 = 0.0;
            for (int j = i + 1; j < n; ++j)
            {
                sum2 += A[i][j] * x[j]; // Используем старые значения
            }

            double new_xi = (b[i] - sum1 - sum2) / A[i][i];
            new_xi = (1 - omega) * x[i] + omega * new_xi;

            diff = std::max(diff, std::abs(new_xi - x[i]));

#pragma omp critical
            x[i] = new_xi; // Обновляем значение переменной
        }

        ++iter;
    } while (diff > TOL && iter < MAX_ITER);

    if (iter == MAX_ITER)
    {
        std::cout << "Метод Верхней Релаксации не сошелся за " << MAX_ITER << " итераций.\n";
    }
    else
    {
        std::cout << "Метод Верхней Релаксации сошелся за " << iter << " итераций.\n";
    }
}

int main()
{
    // Пример СЛАУ
    std::vector<std::vector<double>> A = {
        {5, 2, -7, 14, 0, 0},
        {5, -1, 8, -13, 3, 0},
        {10, 1, -2, 7, -1, 0},
        {15, 3, 15, 9, 7, 0},
        {2, -1, -4, 5, -7, 0}};

    int n = A.size();                               // Размерность матрицы СЛАУ
    std::vector<double> b = {21, 12, 29, 130, -13}; // Вектор правой части
    std::vector<double> x(5, 0.0);                  // Начальные приближения
    double omega = 1.25;                            // Параметр релаксации

    // Прямой метод Гаусса
    gauss_forward_elimination(A, b);
    // Обратная подстановка
    gauss_back_substitution(A, b, x);
    // Вывод решения
    std::cout << "Решение системы:\n";
    for (int i = 0; i < n; ++i)
    {
        std::cout << "x[" << i << "] = " << x[i] << std::endl;
    }

    jacobi_parallel(A, b, x);
    // Вывод результатов метода Якоби
    std::cout << "Решение:\n";
    for (size_t i = 0; i < x.size(); ++i)
    {
        std::cout << "x[" << i << "] = " << x[i] << std::endl;
    }

    gauss_seidel_parallel(A, b, x);
    // Вывод результатов метода Зейделя
    std::cout << "Решение:\n";
    for (size_t i = 0; i < x.size(); ++i)
    {
        std::cout << "x[" << i << "] = " << x[i] << std::endl;
    }

    sor_parallel(A, b, x, omega);
    // Вывод результатов метода Верхней релаксации
    std::cout << "Решение:\n";
    for (size_t i = 0; i < x.size(); ++i)
    {
        std::cout << "x[" << i << "] = " << x[i] << std::endl;
    }

    return 0;
}