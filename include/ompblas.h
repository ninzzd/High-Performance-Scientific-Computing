#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#pragma once
int num_threads;
int printVect(double* x, int n);
int set_num_threads(int n);
int get_num_threads();
int vvdot(double *a, double *b, int n, double *result);
int mvdot(double **a, double *b, int m, int n, double *result);
