#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#define idx(i,j,N) ((j)*(N)+(i))
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
// Size of up, down, left, right = N-2 or M-2 (solution at corner points are always known)
void generateMat(int N, int M, double* down, double* up, double* left, double* right, double x0, double y0, double dx, double dy, double *A, double* b) {
    A = (double*)malloc((N-2) * (M-2) * (N-2) * (M-2) * sizeof(double));
    b = (double*)malloc((N-2) * (M-2) * sizeof(double));
    double coeff_x = 1.0 / (dx * dx);
    double coeff_y = 1.0 / (dy * dy);
    double coeff_center = -2.0 * (coeff_x + coeff_y);
    for (int j = 1; j < M-1; j++) {
        for (int i = 1; i < N-1; i++) {
            int b_idx = idx(i-1, j-1, N-2);
            int a_row = b_idx*(N-2)*(M-2);
            A[a_row + b_idx] = coeff_center;
            if(j==1 && i == 1){
                b[b_idx] = -coeff_x * left[j-1] - coeff_y * down[i-1];
                A[a_row + idx(i-1, j, N-2)] = coeff_x;   // up 
                A[a_row + idx(i, j-1, N-2)] = coeff_y; // right
            }  
            else if(j==1 && i == N-2){
                b[b_idx] = -coeff_x * right[j-1] - coeff_y * down[i-1];
                A[a_row + idx(i-1, j, N-2)] = coeff_x;   // up 
                A[a_row + idx(i-2, j-1, N-2)] = coeff_y; // left
            }  
            else if(j==M-2 && i == 1){
                b[b_idx] = -coeff_x * left[j-1] - coeff_y * up[i-1];
                A[a_row + idx(i-1, j-2, N-2)] = coeff_x;   // down 
                A[a_row + idx(i, j-1, N-2)] = coeff_y; // right
            }  
            else if(j==M-2 && i == N-2){
                b[b_idx] = -coeff_x * right[j-1] - coeff_y * up[i-1];
                A[a_row + idx(i-1, j-2, N-2)] = coeff_x;   // down 
                A[a_row + idx(i-2, j-1, N-2)] = coeff_y; // left
            }  
            else if(i == 1){
                b[b_idx] = -coeff_x * left[j-1];
                A[a_row + idx(i, j-1, N-2)] = coeff_x;   // right
                A[a_row + idx(i-1, j-2, N-2)] = coeff_y; // down
                A[a_row + idx(i-1, j, N-2)] = coeff_y;   // up
            }  
            else if(i == N-2){
                b[b_idx] = -coeff_x * right[j-1];
                A[a_row + idx(i-2, j-1, N-2)] = coeff_x; // left
                A[a_row + idx(i-1, j-2, N-2)] = coeff_y; // down
                A[a_row + idx(i-1, j, N-2)] = coeff_y;   // up
            }  
            else if(j == 1){
                b[b_idx] = -coeff_y * down[i-1];
                A[a_row + idx(i-2, j-1, N-2)] = coeff_x; // left
                A[a_row + idx(i, j-1, N-2)] = coeff_x;   // right
                A[a_row + idx(i-1, j, N-2)] = coeff_y;   // up
            }  
            else if(j == M-2){
                b[b_idx] = -coeff_y * up[i-1];
                A[a_row + idx(i-2, j-1, N-2)] = coeff_x; // left
                A[a_row + idx(i, j-1, N-2)] = coeff_x;   // right
                A[a_row + idx(i-1, j-2, N-2)] = coeff_y; // down
            }  
            else{
                b[b_idx] = 0.0;
                A[a_row + idx(i-2, j-1, N-2)] = coeff_x; // left
                A[a_row + idx(i, j-1, N-2)] = coeff_x;   // right
                A[a_row + idx(i-1, j-2, N-2)] = coeff_y; // down
                A[a_row + idx(i-1, j, N-2)] = coeff_y;   // up
            }
        }
    }
}