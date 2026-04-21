#include "slau.h"
#include <algorithm>
#include <cmath>

void gaussNoPivot(double** A, double* b, int n, double* x) {
    for (int k = 0; k < n - 1; k++) {
        double m = A[k + 1][k] / A[k][k];
        for (int j = 0; j < n; j++) {
            A[k + 1][j] -= A[k][j] * m;
        }
        b[k + 1] -= b[k] * m;
    }

    for (int i = n - 1; i >= 0; i--) {
        double sum = b[i];
        for (int j = i + 1; j < n; j++) {
            sum -= A[i][j] * x[j];
        }
        x[i] = sum / A[i][i];
    }
}

void gaussWithPivot(double** A, double* b, int n, double* x) {
    for (int k = 0; k < n-1; k++) {
        int max_i = k;
        double max = fabs(A[k][k]);

        for (int i = k; i < n; i++) {
            if (fabs(A[i][k]) > max) {
                max = fabs(A[i][k]);
                max_i = i;
            }
        }

        if (max_i != k) {
            for (int i = 0; i < n; i++) {
                std::swap(A[k][i], A[max_i][i]);
            }
            std::swap(b[k], b[max_i]);
        }

        for (int i = k + 1; i < n; i++) {
            double m = A[i][k] / A[k][k];
            for (int j = k; j < n; j++)
                A[i][j] -= m * A[k][j];
            b[i] -= m * b[k];
        }
    }

    for (int i = n - 1; i >= 0; i--) {
        double sum = b[i];
        for (int j = i + 1; j < n; j++) {
            sum -= A[i][j] * x[j];
        }
        x[i] = sum / A[i][i];
    }
}

void LU(double** A, double** L, int* P, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            L[i][j] = (i == j) ? 1.0 : 0.0;
        }
        P[i] = i;
    }

    for (int k = 0; k < n-1; k++) {
        int max_i = k;
        double max = fabs(A[k][k]);
        for (int i = k; i < n; i++) {
            if (fabs(A[i][k]) > max) {
                max = fabs(A[i][k]);
                max_i = i;
            }
        }

        if (max_i != k) {
            for (int j = 0; j < n; j++) {
                std::swap(A[k][j], A[max_i][j]);
                if (j < k) {
                    std::swap(L[k][j], L[max_i][j]);
                }
            }
            std::swap(P[k], P[max_i]);
        }

        for (int i = k + 1; i < n; i++) {
            double m = A[i][k] / A[k][k];
            L[i][k] = m;
            for (int j = k; j < n; j++) {
                A[i][j] -= m * A[k][j];
            }
        }
    }
}

void LowerSubstitution(double** L, double* b, double* y, int n) {
    for (int i = 0; i < n; i++) {
        y[i] = b[i];
        for (int j = 0; j < i; j++) {
            y[i] -= L[i][j] * y[j];
        }
    }
}

void UpperSubstitution(double** U, double* y, double* x, int n) {
    for (int i = n - 1; i >= 0; i--) {
        x[i] = y[i];
        for (int j = i + 1; j < n; j++) {
            x[i] -= U[i][j] * x[j];
        }
        x[i] /= U[i][i];
    }
}