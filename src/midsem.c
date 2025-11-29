#include <iter_solv.h>
int main(){
    int n;
    double **coeff;
    int **jcoeff; 
    double **a, *b;
    printf("Enter n: \n");
    scanf("%d",&n);
    printf("Enter coeffs: \n");

    coeff = (double**)malloc(n*sizeof(double*));
    jcoeff = (int**)malloc(n*sizeof(int*));
    a = (double**)malloc(n*sizeof(double*));
    b = (double*)malloc(n*sizeof(double));
    for(int i = 0;i < n;i++){
        coeff[i] = (double*)malloc(3*sizeof(double));
        for(int j = 0;j < 3;j++){
            scanf("%lf",&coeff[i][j]);
        }
    }
    printf("Enter jcoeffs: \n");
    for(int i = 0;i < n;i++){
        jcoeff[i] = (int*)malloc(3*sizeof(int));
        for(int j = 0;j < 3;j++){
            scanf("%d",&jcoeff[i][j]);
        }
    }
    // printf("Coeffs:\n");
    // for(int i = 0;i < n;i++){
    //     for(int j = 0;j < 3;j++){
    //         printf("%lf ",coeff[i][j]);
    //     }
    //     printf("\n");

    // }
    // printf("Jcoeffs:\n");
    // for(int i = 0;i < n;i++){
    //     for(int j = 0;j < 3;j++){
    //         printf("%d ",jcoeff[i][j]);
    //     }
    //     printf("\n");

    // }
    for(int i = 0;i < n;i++){
        a[i] = (double*) malloc(n*sizeof(double));
        for(int j = 0;j < n;j++){
            a[i][j] = 0.0;
        }
    }
    for(int i = 0;i < n;i++){
        for(int j = 0;j < 3;j++){
            if(coeff[i][j] != 0.0){
                int col = jcoeff[i][j] - 1;
                a[i][col] = coeff[i][j];
            }
        }
    }
    // printf("Matrix A obtained: \n");
    // for(int i = 0;i < n;i++){
    //     for(int j = 0;j < n;j++){
    //         printf("%lf ",a[i][j]);
    //     }
    //     printf("\n");
    // }   
    printf("Enter b:\n");
    for(int i = 0;i < n;i++){
        scanf("%lf",&b[i]);
    }
    double *res1;
    res1 = (double*)malloc(n*sizeof(double));
    double *x0 = (double*)malloc(n*sizeof(double));
    for(int i = 0;i < n;i++){
        x0[i] = 0.0;
    }
    double w;
    printf("Enter the SOR factor: \n");
    scanf("%lf",&w);
    int sor_iter = sorSolver(a,b,x0,n,0.000001,w,1,res1);
    // int sor_2 = sorSolver(a,b,x0,n,0.000001,1.6,1,res2);
    // printf("No. of iterations for w = 1.6: %d\n",sor_2);
    // printVect(res2,n);
    return 0;
}