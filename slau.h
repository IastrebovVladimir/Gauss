#ifndef SLAU_H
#define SLAU_H

void gaussNoPivot(double** A, double* b, int n, double* x);
void gaussWithPivot(double** A, double* b, int n, double* x);
void LU(double** A, double** L, int* P, int n);
void LowerSubstitution(double** L, double* b, double* y, int n);
void UpperSubstitution(double** U, double* y, double* x, int n);

#endif
