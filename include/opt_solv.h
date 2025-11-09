#include "cblas.h"
#include "lapacke.h"
#include "stdio.h"
#include "stdlib.h"
#include "math.h"
#include "iter_solv.h"
#include "float.h"
#include "omp.h"
#pragma once
int dot(double *a, double *b, double *res, int n);
double l2nrm(double *a, int n);
int steepestDescent(double **a, double *b, double* x0, int n, double eta, int mode, double* result);
int minimalResidual(double **a, double *b, double* x0, int n, double eta, int mode, double* result);
int conjugateGradient(double **a, double *b, double *x0, int n, double eta, int mode, double *result);
int bicgstab(double **a, double *b, double *x0, int n, double eta, int mode, double *result);
int bicgstabOMP(double **a, double *b, double *x0, int n, double eta, int mode, double *result);