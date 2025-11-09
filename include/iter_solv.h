#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <lapacke.h>
#include <cblas.h>
#include <omp.h>
#include <ompblas.h>
#pragma once
double sorOpt(double **a, int n);
long long int jacobiSolver(double** a, double* b, double* x_0, int n, double eta, int mode, double* result);
long long int gaussSeidelSolver(double** a, double* b, double* x_0, int n, double eta, int mode, double* result);
long long int sorSolver(double** a, double* b, double* x_0, int n, double eta, double w, int mode, double* result);
void printIter(long long int iter);
// void printVect(double* v, int n);
void printMat(double** a, double* b, double* x, int n);
void fGetMat(FILE** kmat, FILE** fvec, FILE** kinfo);
void getMat(double*** a, double** b,int* n, FILE* kmat, FILE* fvec, FILE* kinfo);
