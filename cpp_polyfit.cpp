#include <iostream>
#include <math.h>

using namespace std;


// Function prototype
int polyfit(const double* const dependentValues,
            const double* const independentValues,
            unsigned int        countOfElements,
            unsigned int        order,
            double*             coefficients);

// Constant values
const double Power[6] = {120.0, 80.0, 60.0, 20.0, 5.0, 3.0};
const double Distance[6] = {0.0, 25.0, 50.0, 100.0, 200.0, 250.0};


// Main function
int main()
{
  cout << "This is the universal polyfit test program.\n\n";

  int N = 3; // степень экстраполяции
  double coefs[N+1];
  polyfit(Distance, Power, 6, N, coefs);
  cout << "Polynom coefs: " << coefs[0] << "\t" << coefs[1] << "\t" << coefs[2] << "\t" << coefs[3] << endl;
  return 0;
}


// Функция для экстрополяции данных полиномом заданной степени
int polyfit(const double* const dependentValues,   // зависимые величины
            const double* const independentValues, // независимые величины
            unsigned int        countOfElements,   // количество элементов в массиве данных
            unsigned int        order,             // степень полинома для экстраполяции
            double*             coefficients)      // выходной массив коэффициентов
{
    // Объявление переменных
    // ----------------------------------
    enum {maxOrder = 5};

    double B[maxOrder+1] = {0.0f};
    double P[((maxOrder+1) * 2)+1] = {0.0f};
    double A[(maxOrder + 1)*2*(maxOrder + 1)] = {0.0f};

    double x, y, powx;

    unsigned int ii, jj, kk;

    // Проверка начальных условий
    // ----------------------------------

    // Необходимое условие: countOfElements > (order+1)
    if (countOfElements <= order)
        return -1;

    // Данная функция эквивалентна функции polyfit из
    // среды Matlab при степени эстраполяции <= 5.
    // При более высоких степенях появляются расхождения в решении.
    if (order > maxOrder)
        return -1;

    // Начало программы:
    // ----------------------------------

    // Нахождение вектора-столбца
    for (ii = 0; ii < countOfElements; ii++)
    {
        x    = dependentValues[ii];
        y    = independentValues[ii];
        powx = 1;

        for (jj = 0; jj < (order + 1); jj++)
        {
            B[jj] = B[jj] + (y * powx);
            powx  = powx * x;
        }
    }

    // Создание массива PowX
    P[0] = countOfElements;

    // Вычисление сумм степеней элементов массива Х
    for (ii = 0; ii < countOfElements; ii++)
    {
        x    = dependentValues[ii];
        powx = dependentValues[ii];

        for (jj = 1; jj < ((2 * (order + 1)) + 1); jj++)
        {
            P[jj] = P[jj] + powx;
            powx  = powx * x;
        }
    }

    // Заполнение матрицы А (reduction matrix)
    for (ii = 0; ii < (order + 1); ii++)
    {
        for (jj = 0; jj < (order + 1); jj++)
        {
            A[(ii * (2 * (order + 1))) + jj] = P[ii+jj];
        }

        A[(ii*(2 * (order + 1))) + (ii + (order + 1))] = 1;
    }

    // Перенос единичной части матрицы А на левую сторону
    // (нахождение обратной матрицы для левой стороны матрицы А)
    for (ii = 0; ii < (order + 1); ii++)
    {
        x = A[(ii * (2 * (order + 1))) + ii];
        if (x != 0)
        {
            for (kk = 0; kk < (2 * (order + 1)); kk++)
            {
                A[(ii * (2 * (order + 1))) + kk] =
                    A[(ii * (2 * (order + 1))) + kk] / x;
            }

            for (jj = 0; jj < (order + 1); jj++)
            {
                if ((jj - ii) != 0)
                {
                    y = A[(jj * (2 * (order + 1))) + ii];
                    for (kk = 0; kk < (2 * (order + 1)); kk++)
                    {
                        A[(jj * (2 * (order + 1))) + kk] =
                            A[(jj * (2 * (order + 1))) + kk] -
                            y * A[(ii * (2 * (order + 1))) + kk];
                    }
                }
            }
        }
        else
        {
            // Для вырожденной матрицы решения нет
            return -1;
        }
    }

    // Расчёт коэффициентов полинома
    for (ii = 0; ii < (order + 1); ii++)
    {
        for (jj = 0; jj < (order + 1); jj++)
        {
            x = 0;
            for (kk = 0; kk < (order + 1); kk++)
            {
                x = x + (A[(ii * (2 * (order + 1))) + (kk + (order + 1))] *
                    B[kk]);
            }
            coefficients[ii] = x;
        }
    }

    return 0;
}
