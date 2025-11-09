#include <omp.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>
#pragma once
int n_threads;
int printVect(double* x, int n);
int printCPUTime(struct timeval start, struct timeval end);
int set_num_threads(int n);
int get_num_threads();
int vvdot(double *a, double *b, int n, double *result);
int mvdot(double **a, double *b, int m, int n, double *result);
