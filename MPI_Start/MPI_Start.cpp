#include <iostream>
#include <omp.h>
#include <cstdlib>
#include <chrono>
using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
using std::chrono::milliseconds;
using namespace std;
int main(int n, int t) {
    n = 1000;
    t = 15;
    setlocale(LC_ALL, "RUS");
    double** matrix = new double* [n];
    for (int i = 0; i < n; i++)
        matrix[i] = new double[n + 1];
    
    // Инициализация матрицы
    for (int i = 0; i < n; i++) {
        double b;
        for (int j = 0; j < n + 1; j++) {
            matrix[i][j] = rand() % 10;
            /*if (j == n) {
                b = matrix[i][j];
                
            }
            else if (j + 1 == n) {
                cout << matrix[i][j] << "x" << j + 1;
            }
            else {
                cout << matrix[i][j] << "x" << j + 1;
                cout << " + ";
            }*/
        }
        // cout << " = " << b << endl;
    }

    omp_set_num_threads(t);
    double timein = omp_get_wtime();

    // Прямой ход метода Гаусса
    for (int i = 0; i < n; i++) {
        // Параллельный блок для обновления строк ниже текущей
#pragma omp parallel for
        for (int j = i + 1; j < n; j++) {
            double factor = static_cast<double>(matrix[j][i]) / matrix[i][i];
            for (int k = i; k <= n; k++)
                matrix[j][k] -= factor * matrix[i][k];
        }
    }

    // Обратный ход
    double* xx = new double[n];
    xx[n - 1] = matrix[n - 1][n] / matrix[n - 1][n - 1];
    for (int i = n - 2; i >= 0; i--) {
        xx[i] = matrix[i][n];
        for (int j = i + 1; j < n; j++)
            xx[i] -= matrix[i][j] * xx[j];
        xx[i] /= matrix[i][i]; // делим на коэффициент перед x[i]
    }

    double timeout = omp_get_wtime();
    double dt = timeout - timein;

    /*cout << "\nРешение системы уравнений:";
    for (int i = 0; i < n; i++)
        cout << "\nx" << i + 1 << " = " << xx[i];*/

    std::cout << endl << "Время вычислений : " << dt*1000 << " милисекунд" << std::endl;

    /*auto end = high_resolution_clock::now();

    
    duration<double, std::milli> ms_double = end - start;
    std::cout << endl << ms_double.count() << "ms\n";*/

    delete[] xx;
    for (int i = 0; i < n; i++)
        delete[] matrix[i];
    delete[] matrix;
}