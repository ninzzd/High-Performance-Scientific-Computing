#include <ompblas.h>
int n_threads = 1;
double randd(double min, double max){
    return min + ((double)rand()/(double)RAND_MAX)*(max-min);
}
int printVect(double* x, int n){
    if(n < 0)
        return -1;
    printf("[");
    for(int i = 0;i < n;i++){
        if(i == 0)
        printf("%lf",x[i]);
        else
            printf(" %lf",x[i]);
    }
    printf("]\n");
    return 0;
}
int getCPUTime(struct timeval start, struct timeval end, int mode){
    time_t s, us, ms;
    s = end.tv_sec - start.tv_sec;
    us = end.tv_usec - start.tv_usec;
    if(us < 0){
        us += 1000000;
        s -= 1;
    }
    ms = us/1000;
    if (mode == 1)
        printf("CPU Time = %lds %ldms %ldus\n",s,ms,us%1000);
    return (int)(s*1e6 + us);
}
int set_num_threads(int n){
    if(n < 1)
        return -1;
    omp_set_num_threads(n);
    n_threads = n;
    return 0;
}
int get_num_threads(){
    return n_threads;
}
int vvdot(double *a, double *b, int n, double *result){
    *result = 0;    
    for(int i = 0;i < n;i++){
        *result += a[i]*b[i];
    }
    return 0;
}
int vvdot_omp(double *a, double *b, int n, double *result){
    *result = 0;
    // printf("Number of threads = %d\n",n_threads);
    double *partial = (double*)malloc(n*sizeof(double));
    // Parallelized multiplication only
    int i;
    #pragma omp parallel for num_threads(n_threads) private(i) default(shared)
    for(i = 0;i < n;i++){
        partial[i] = a[i]*b[i];
        // id = omp_get_thread_num();
        // printf("\ti = %d, id = %d\n",i,id);
    }
    
    for(i = 0;i < n;i++){
        *result += partial[i];
    }
    free(partial);
    return 0;
}
int mvdot(double **a, double *b, int m, int n, double *res){
    // result should be allocated memory in the caller scope
    int i,j;
    for(i = 0;i < m;i++){
        res[i] = 0.0f;
        for(j = 0;j < n;j++){
            res[i] += a[i][j]*b[j];
        }
    }
    return 0;
}
int mvdot_omp(double **a, double *b, int m, int n, double *res){
    // result should be allocated memory in the caller scope
    int i,j;
    double temp;
    #pragma omp parallel for num_threads(n_threads) private(i,j,temp) default(shared)
    for(i = 0;i < m;i++){
        temp = 0.0;
        for(j = 0;j < n;j++){
            temp += a[i][j]*b[j];
        }
        res[i] = temp; // Race condition does occur but all threads with same value of 'i' are trying to force the same value of 'temp' to res[i]
    }
    return 0;
}
int dnrm2(double *a, int n, double *res){
    *res = 0.0;
    for(int i = 0;i < n;i++){
        *res += a[i]*a[i];
    }
    *res = sqrt(*res);
    return 0;
}
// int dnrm2_omp(double *a, int n, double *res){
//     *res = 0.0;
//     for(int i = 0;i < n;i++){
//         res += a[i]*a[i];
//     }
//     res = sqrt(res);
//     return 0;
// }