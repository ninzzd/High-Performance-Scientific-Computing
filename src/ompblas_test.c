#include <ompblas.h>
#include <cblas.h>
#include <math.h>
#include "sys/time.h"
#include "time.h"
#include "unistd.h"
double randd(double min, double max){
    return min + ((double)rand()/(double)RAND_MAX)*(max-min);
}
int main(){
    FILE *bmark; // Benchmark
    bmark = fopen("ompblas_bench.csv","w");
    struct timeval start,end;
    int n,n_max;
    n_max = 1000;
    // printf("Enter n = ");
    // scanf("%d",&n);
    // printf("\n");
    double min,max;
    min = -1.0;
    max = 1.0;
    double *a, *b, res1, res2;
    for(n = 2;n <= n_max;n++){
        printf("n = %d\n",n);
        a = (double*)malloc(n*sizeof(double)); 
        b = (double*)malloc(n*sizeof(double)); 
        for(int i = 0;i < n;i++){
            a[i] = randd(min,max);
            b[i] = randd(min,max);
        }
        // printVect(a,n);
        // printVect(b,n);

        gettimeofday(&start,NULL);
        vvdot(a,b,n,&res1);
        gettimeofday(&end,NULL);

        printf("Result (Using ompblas) = %lf\n",res1);
        int ompblas = printCPUTime(start,end);

        gettimeofday(&start,NULL);
        res2 = cblas_ddot(n,a,1,b,1);
        gettimeofday(&end,NULL);
        printf("Result (Using cblas) = %lf\n",res2);
        int cblas = printCPUTime(start,end);
        printf("Is correct = %d\n",fabs(res1 - res2) <= 1e-6);
        fprintf(bmark,"%d,%d,%d\n",n,ompblas,cblas);
    }
    return 0;
}