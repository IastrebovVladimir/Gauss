#include <iostream>
#include <iomanip>
#include <random>
#include <chrono>
#include "slau.h"

int main() {
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "n\tGauss\tLU(b1)\tLU(b2)\n";

    int sizes[] = {100, 200, 500, 1000};

    std::mt19937 gen(std::random_device{}());
    std::uniform_real_distribution<double> dis(-1.0, 1.0);

    for (int n : sizes) {
        double** A = new double*[n];
        double** L = new double*[n];
        double** A_g = new double*[n];
        double** A_lu = new double*[n];

        for (int i = 0; i < n; i++) {
            A[i] = new double[n];
            L[i] = new double[n];
            A_g[i] = new double[n];
            A_lu[i] = new double[n];
        }

        double* b1 = new double[n];
        double* b2 = new double[n];
        double* b_g = new double[n];
        double* b1_lu = new double[n];
        double* b2_lu = new double[n];

        double* x_gauss = new double[n];
        double* x_lu1 = new double[n];
        double* x_lu2 = new double[n];
        double* y = new double[n];
        int* P = new int[n];

        // Генерация матрицы A
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                A[i][j] = dis(gen);
                if (i == j) {
                    A[i][j] += 2.0;
                }
            }
        }

        // Генерация правых частей
        for (int i = 0; i < n; i++) {
            b1[i] = dis(gen);
            b2[i] = dis(gen);
        }

        // Копия для метода Гаусса
        for (int i = 0; i < n; i++) {
            b_g[i] = b1[i];
            for (int j = 0; j < n; j++) {
                A_g[i][j] = A[i][j];
            }
        }

        auto t1 = std::chrono::high_resolution_clock::now();
        gaussWithPivot(A_g, b_g, n, x_gauss);
        auto t2 = std::chrono::high_resolution_clock::now();
        double gauss_time = std::chrono::duration<double, std::milli>(t2 - t1).count();

        // Копия для LU-разложения
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                A_lu[i][j] = A[i][j];
                L[i][j] = (i == j) ? 1.0 : 0.0;
            }
        }

        auto t3 = std::chrono::high_resolution_clock::now();
        LU(A_lu, L, P, n);
        auto t4 = std::chrono::high_resolution_clock::now();
        double lu_time = std::chrono::duration<double, std::milli>(t4 - t3).count();

        // Решение для b1 через уже готовое LU
        for (int i = 0; i < n; i++) {
            b1_lu[i] = b1[P[i]];
        }

        auto t5 = std::chrono::high_resolution_clock::now();
        LowerSubstitution(L, b1_lu, y, n);
        UpperSubstitution(A_lu, y, x_lu1, n);
        auto t6 = std::chrono::high_resolution_clock::now();
        double solve1_time = std::chrono::duration<double, std::milli>(t6 - t5).count();

        // Решение для b2 через то же LU-разложение
        for (int i = 0; i < n; i++) {
            b2_lu[i] = b2[P[i]];
        }

        auto t7 = std::chrono::high_resolution_clock::now();
        LowerSubstitution(L, b2_lu, y, n);
        UpperSubstitution(A_lu, y, x_lu2, n);
        auto t8 = std::chrono::high_resolution_clock::now();
        double solve2_time = std::chrono::duration<double, std::milli>(t8 - t7).count();

        std::cout << n << "\t"
                  << gauss_time << "\t"
                  << (lu_time + solve1_time) << "\t"
                  << solve2_time << "ms\n";

        for (int i = 0; i < n; i++) {
            delete[] A[i];
            delete[] L[i];
            delete[] A_g[i];
            delete[] A_lu[i];
        }

        delete[] A;
        delete[] L;
        delete[] A_g;
        delete[] A_lu;

        delete[] b1;
        delete[] b2;
        delete[] b_g;
        delete[] b1_lu;
        delete[] b2_lu;

        delete[] x_gauss;
        delete[] x_lu1;
        delete[] x_lu2;
        delete[] y;
        delete[] P;
    }

    return 0;
}