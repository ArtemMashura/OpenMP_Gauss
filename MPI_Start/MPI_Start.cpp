//#include <iostream>
//#include <time.h>
//#include <omp.h>
//#include <vector>
//#include <chrono>
//using std::chrono::high_resolution_clock;
//using std::chrono::duration_cast;
//using std::chrono::duration;
//using std::chrono::milliseconds;
//using namespace std;
//
//bool search_reverse_matrix(vector <vector<double>>& matrix)
//{
//    int size = matrix.size();
//    vector <vector<double>> E(size, vector<double>(size));
//
//    //Заполнение единичной матрицы
//    for (int i = 0; i < size; i++)
//    {
//        for (int j = 0; j < size; j++)
//        {
//            if (i == j) E[i][j] = 1.0;
//            else E[i][j] = 0.0;
//        }
//    }
//
//    //Трансформация исходной матрицы в верхнетреугольную
//    for (int k = 0; k < size; k++)
//    {
//        if (abs(matrix[k][k]) < 1e-8)
//        {
//            bool changed = false;
//
//            for (int i = k + 1; i < size; i++)
//            {
//                if (abs(matrix[i][k]) > 1e-8)
//                {
//                    swap(matrix[k], matrix[i]);
//                    swap(E[k], E[i]);
//                    changed = true;
//                    break;
//                }
//            }
//
//            if (!changed)
//                return false;
//        }
//
//        double div = matrix[k][k];
//
//#pragma omp parallel
//        {
//#pragma omp for
//            for (int j = 0; j < size; j++)
//            {
//                matrix[k][j] /= div;
//                E[k][j] /= div;
//            }
//        }
//
//#pragma omp parallel
//        {
//#pragma omp for
//            for (int i = k + 1; i < size; i++)
//            {
//                double multi = matrix[i][k];
//
//
//                for (int j = 0; j < size; j++)
//                {
//                    matrix[i][j] -= multi * matrix[k][j];
//                    E[i][j] -= multi * E[k][j];
//                }
//            }
//        }
//    }
//
//    //Формирование единичной матрицы из исходной
//    //и обратной из единичной
//    for (int k = size - 1; k > 0; k--)
//    {
//#pragma omp parallel
//        {
//#pragma omp for
//            for (int i = k - 1; i > -1; i--)
//            {
//                double multi = matrix[i][k];
//
//                for (int j = 0; j < size; j++)
//                {
//                    matrix[i][j] -= multi * matrix[k][j];
//                    E[i][j] -= multi * E[k][j];
//                }
//            }
//        }
//    }
//    matrix = E;
//    return true;
//}
//
//int t = 10;
//int main()
//{
//    setlocale(LC_ALL, "RUS");
//    int equations_amount;
//    cout << "Введите количество уравнений: ";
//    cin >> equations_amount;
//
//    vector<vector<double>> matrix(equations_amount, vector<double>(equations_amount));
//    vector<double> B(equations_amount);
//
//    // Заполняем матрицу коэффициентов и B
//    for (int i = 0; i < equations_amount; i++)
//    {
//        for (int j = 0; j < equations_amount; j++) {
//            matrix[i][j] = rand() % 10;
//            /*cout << matrix[i][j] << "x" << j + 1;
//            if (j + 1 < equations_amount) {
//                cout << " + ";
//            }*/
//        }
//
//        B[i] = rand() % 10;
//        // cout << " = " << B[i] << endl;
//    }
//
//    // double t = clock();
//    omp_set_num_threads(t);
//    auto start = high_resolution_clock::now();
//    // Вычисление обратной матрицы
//    if (!search_reverse_matrix(matrix))
//    {
//        cout << "\nНевозможно найти обратную матрицу!\n";
//        exit(1);
//    }
//
//    // Матрица-столбец неизвестных X и вычисление окончательного результата
//    vector<double> X(equations_amount);
//#pragma omp parallel
//    {
//#pragma omp for
//        for (int i = 0; i < equations_amount; i++)
//        {
//            X[i] = 0;
//            for (int j = 0; j < equations_amount; j++)
//                X[i] += matrix[i][j] * B[j];
//        }
//    }
//
//    // Вывод окончательного результата
//   /* cout << "\nРешение системы уравнений:";
//    for (int i = 0; i < equations_amount; i++)
//        cout << "\nx" << i + 1 << " = " << X[i];*/
//
//    // t = (clock() - t) / 1000;
//    // cout << "\n\nВремя, затраченное на вычисление: " << t << "с.\n";
//    auto end = high_resolution_clock::now();
//    duration<double, std::milli> ms_double = end - start;
//    std::cout << endl << ms_double.count() << "ms\n";
//    return 0;
//}









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