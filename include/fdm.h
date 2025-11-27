#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#pragma once
void generateFDM();
void generateMat(int N, int M, double *A, double* b, double* down, double* up, double* left, double* right, double x0, double y0, double dx, double dy);