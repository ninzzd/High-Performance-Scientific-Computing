#include <ompblas.h>
#include <cblas.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>
double randd(double min, double max){
    return min + ((double)rand()/(double)RAND_MAX)*(max-min);
}
int main(){
    FILE *bmark; // Benchmark
    bmark = fopen("ompblas_bench.csv","w");
    struct timeval start,end;
    int n,n_max;
    n_max = 1e4;
    // n = 5;
    // printf("Enter n = ");
    // scanf("%d",&n);
    // printf("\n");
    set_num_threads(8);
    fprintf(bmark,"n_threads,%d,,\n",get_num_threads());
    fprintf(bmark,"n,omp,cblas,for\n");
    double min,max;
    min = -1.0;
    max = 1.0;
    double *a, *b, res1, res2, res3;
    // a = (double*)malloc(n*sizeof(double)); 
    // b = (double*)malloc(n*sizeof(double)); 
    // for(int i = 0;i < n;i++){
    //     a[i] = randd(min,max);
    //     b[i] = randd(min,max);
    // }
    // vvdot(a,b,n,&res1);
    for(n = 2;n <= n_max;n++){
        res1 = res2 = res3 = 0.0;
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

        // printf("Result (Using ompblas) = %lf\n",res1);
        int ompblas = printCPUTime(start,end);

        gettimeofday(&start,NULL);
        // res2 = cblas_ddot(n,a,1,b,1);
        for(int i = 0;i < n;i++){
            res2 += a[i]*b[i];
        }
        gettimeofday(&end,NULL);
        // printf("Result (Using cblas) = %lf\n",res2);
        int iter = printCPUTime(start,end);
        // printf("Is correct = %d\n",fabs(res1 - res2) <= 1e-6);

        gettimeofday(&start,NULL);
        res3 = cblas_ddot(n,a,1,b,1);
        gettimeofday(&end,NULL);
        int cblas = printCPUTime(start,end);
        fprintf(bmark,"%d,%d,%d,%d\n",n,ompblas,cblas,iter);
        free(a);
        free(b);
    }
    return 0;
}