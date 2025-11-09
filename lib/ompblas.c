#include <ompblas.h>
int num_threads = 1;
int printVect(double* x, int n){
    if(n < 0)
        return -1;
    printf("[");
    for(int i = 0;i < n;i++){
        printf("%lf ",x[i]);
    }
    printf("]\n");
    return 0;
}
int set_num_threads(int n){
    if(n < 1)
        return -1;
    omp_set_num_threads(n);
    num_threads = n;
    return 0;
}
int get_num_threads(){
    return num_threads;
}
int vvdot(double *a, double *b, int n, double *result){
    *result = 0;
    double *partial = (double*)malloc(n*sizeof(double));
    // Parallelized multiplication only
    #pragma omp parallel for default(shared) private(i)
    {
        for(int i = 0;i < n;i++){
            partial[i] = a[i]*b[i];
        }
    }
    for(int i = 0;i < n;i++){
        *result += partial[i];
    }
    return 0;
}
int mvdot(double **a, double *b, int m, int n, double *result){
    // result should be allocated memory in the caller scope
    return 0;
}