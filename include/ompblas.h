#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>
#pragma once
int n_threads;
int printVect(double* x, int n);
int printCPUTime(struct timeval start, struct timeval end);
int set_num_threads(int n);
int get_num_threads();
double randd(double min, double max);
int vvdot(double *a, double *b, int n, double *result);
int vvdot_omp(double *a, double *b, int n, double *result);
int mvdot(double **a, double *b, int m, int n, double *result);
int mvdot_omp(double **a, double *b, int m, int n, double *result);
int dnrm2(double *a, int n, double *res);
