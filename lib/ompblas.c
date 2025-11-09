#include <ompblas.h>
int n_threads = 1;
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
int printCPUTime(struct timeval start, struct timeval end){
    time_t s, us, ms;
    s = end.tv_sec - start.tv_sec;
    us = end.tv_usec - start.tv_usec;
    if(us < 0){
        us += 1000000;
        s -= 1;
    }
    ms = us/1000;
    printf("CPU Time = %lds %ldms %ldus\n",s,ms,us%1000);
    return (int)(s*1e6 + ms*1e3 + us);
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
int mvdot(double **a, double *b, int m, int n, double *result){
    // result should be allocated memory in the caller scope
    return 0;
}