#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#pragma once
void generateFDM();
void generateMat(int N, int M, double* down, double* up, double* left, double* right, double dx, double dy, double *A, double* b);