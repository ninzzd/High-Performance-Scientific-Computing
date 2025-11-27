#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
void generateFDM(){
    double d; // Discretization error in both x and y coordinates
    double u0y, u1y, ux0, ux1;
    int n; // No. of points per dimension
    int n_; // Total number of interior points
    double *b;
    double **a;
    char cwd[512] = __FILE__;
    char kmat_str[512];
    char kinfo_str[512];
    char fvec_str[512];
    int null_index;
    // To find the parent directory of the source file
    for(int i = 511;i >= 0;i--){
        if(cwd[i] == '\0'){
            null_index = i;
            break;
        }
    }
    #ifdef _WIN32
    for(int j = null_index-1;j >= 0;j--){
        if(cwd[j] == '\\'){
            break;
        }
        else{
            cwd[j] = '\0';
        }
    }
    #elif __linux__
    for(int j = null_index-1;j >= 0;j--){
        if(cwd[j] == '/'){
            break;
        }
        else{
            cwd[j] = '\0';
        }
    }
    #endif
    FILE *kmat, *fvec, *kinfo;
    printf("----- Laplace Equation Solver -----\n");
    printf("Finite Difference Method: Discretization\n");
    printf("Consider the surface to be a square bounded by the points: (0,0), (0,1), (1,0), (1,1)\n");
    printf("Enter the required approximate discretization error: dx = dy = ");
    scanf("%lf",&d);
    n = 1 + (int)floor(1.0/d); // No. of evenly spaced grid points = No. of divisions + 1
    d = 1.0/(double)(n-1);
    n_ = (n-2)*(n-2);
    b = (double*)malloc(n_*sizeof(double));
    a = (double**) malloc(n_*n_*sizeof(double));
    for(int i = 0;i < n_;i++){
        a[i] = (double*)malloc(n_*sizeof(double));
        for(int j = 0;j < n_;j++){
            a[i][j] = 0.0;
        }
    }
    printf("Grid size: %dx%d\n",n,n);
    printf("Actual discretization error: %lf\n",d);
    printf("Enter the dirichlet boundary conditions:\n");
    printf("For all 0 < y < 1, u(0,y) = ");
    scanf("%lf",&u0y);
    printf("\nFor all 0 < y < 1, u(1,y) = ");
    scanf("%lf",&u1y);
    printf("\nFor all 0 <= x <= 1, u(x,0) = ");
    scanf("%lf",&ux0);
    printf("\nFor all 0 <= x <= 1, u(x,1) = ");
    scanf("%lf",&ux1);
    for(int j = 1;j <= n-2;j++){
        for(int i = 1;i <= n-2;i++){
            int idx_ = (j-1)*(n-2) + (i-1); // index corresponsding to b vector
            if(i == 1 && j == 1){
                b[idx_] = -ux0 - u0y;
                a[idx_][idx_+1] = 1.0; // x + d
                a[idx_][idx_+(n-2)] = 1.0; // y + d           
            }
            else if(i == n-2 && j == n-1){
                b[idx_] = -ux1 - u1y;
                a[idx_][idx_-1] = 1.0; // x - d
                a[idx_][idx_-(n-2)] = 1.0; // y - d             
            }
            else if(i == 1 && j == n-2){
                b[idx_] = -ux1 - u0y;   
                a[idx_][idx_+1] = 1.0; // x + d
                a[idx_][idx_-(n-2)] = 1.0; // y - d           
            }
            else if(i == n-2 && j == 1){
                b[idx_] = -ux0 - u1y;
                a[idx_][idx_-1] = 1.0; // x - d
                a[idx_][idx_+(n-2)] = 1.0; // y + d   
            }
            else if(i == 1){
                b[idx_] = -u0y;
                a[idx_][idx_+1] = 1.0; // x + d
                a[idx_][idx_-(n-2)] = 1.0; // y - d
                a[idx_][idx_+(n-2)] = 1.0; // y + d

            }
            else if(i == n-2){
                b[idx_] = -u1y;
                a[idx_][idx_-1] = 1.0; // x - d
                a[idx_][idx_-(n-2)] = 1.0; // y - d   
                a[idx_][idx_+(n-2)] = 1.0; // y + d   
            }
            else if(j == 1){
                b[idx_] = -ux0;
                a[idx_][idx_-1] = 1.0; // x - d
                a[idx_][idx_+1] = 1.0; // x + d   
                a[idx_][idx_+(n-2)] = 1.0; // y + d   
            }
            else if(j == n-2){
                b[idx_] = -ux1;
                a[idx_][idx_-1] = 1.0; // x - d
                a[idx_][idx_+1] = 1.0; // x + d   
                a[idx_][idx_-(n-2)] = 1.0; // y - d 
            }
            else{
                b[idx_] = 0.0;
                a[idx_][idx_-1] = 1.0; // x - d
                a[idx_][idx_+1] = 1.0; // x + d   
                a[idx_][idx_-(n-2)] = 1.0; // y - d 
                a[idx_][idx_+(n-2)] = 1.0; // y + d
            }
            a[idx_][idx_] = -4.0;
        }
    }
    strcpy(kmat_str,cwd);
    strcpy(kinfo_str,cwd);
    strcpy(fvec_str,cwd);
    strcat(kmat_str,"kmat.txt");
    strcat(kinfo_str,"kinfo.txt");
    strcat(fvec_str,"Fvec.txt");

    kinfo = fopen(kinfo_str,"w");
    fprintf(kinfo,"%d\n",n_);
    fclose(kinfo);
    fvec = fopen(fvec_str,"w");
    for(int i = 0;i < n_;i++){
        fprintf(fvec,"%lf\n",b[i]);
    }
    fclose(fvec);
    kmat = fopen(kmat_str,"w");
    for(int i = 0;i < n_;i++){
        for(int j = 0;j < n_;j++){
            fprintf(kmat,"%lf\n",a[i][j]);
        }
    }
    fclose(kmat);
    // free(b);
    // for(int i = 0;i < n_;i++){
    //     free(a[i]);
    // }
    // free(a);
}
void generateMat(int N, int M, double *A, double* b, double* down, double* up, double* left, double* right, double x0, double y0, double dx, double dy) {
    A = (double*)malloc((N-2) * (M-2) * sizeof(double));
    for (int j = 0; j < M; j++) {
        for (int i = 0; i < N; i++) {
            int idx = j * N + i;
            
        }
    }
}