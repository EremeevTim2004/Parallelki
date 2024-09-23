#include "includes.h"

void bit_reversal(std::vector<int> &data)
{
    int n = data.size();
    int log_n = log2(n);

#pragma omp parallel for
    for (int i = 0; i < n; ++i)
    {
        int reversed_index = 0;
        for (int j = 0; j < log_n; ++j)
        {
            if (i & (1 << j))
            {
                reversed_index |= (1 << (log_n - 1 - j));
            }
        }
        // Запись значения в новом индексе
        if (i < reversed_index)
        {
            std::swap(data[i], data[reversed_index]);
        }
    }
}

int main()
{
    int n = 16; // Размер массива
    std::vector<int> data(n);

    // Заполнение массива
    for (int i = 0; i < n; ++i)
    {
        data[i] = i;
    }

    bit_reversal(data);

    // Вывод результата
    std::cout << "Бит-реверсированный массив: ";
    for (int val : data)
    {
        std::cout << val << " ";
    }
    std::cout << std::endl;

    return 0;
}