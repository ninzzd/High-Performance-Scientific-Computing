#include <ompblas.h>
#include <cblas.h>
#include <math.h>
#include <sys/time.h>
#include <time.h>
#include <unistd.h>
int main(){
    FILE *bmark; // Benchmark
    bmark = fopen("ompblas_mvdot_bench.csv","w");
    struct timeval start,end;
    int n,n_max;
    n_max = 3e3;
    set_num_threads(8);
    fprintf(bmark,"n_threads,%d,\n",get_num_threads());
    fprintf(bmark,"n,omp,no_omp\n");
    double min,max;
    min = -1.0;
    max = 1.0;
    double **a, *b, *res1, *res2,*err;
    for(n = 2;n <= n_max;n++){
        a = (double**)malloc(n*sizeof(double*));
        b = (double*)malloc(n*sizeof(double));
        res1 = (double*)malloc(n*sizeof(double));
        res2 = (double*)malloc(n*sizeof(double));
        err = (double*)malloc(n*sizeof(double));
        for(int i = 0;i < n;i++){
            a[i] = (double*)malloc(n*sizeof(double));
            for(int j = 0;j < n;j++){
                a[i][j] = randd(min,max);
            }
            b[i] = randd(min,max);
        }
        gettimeofday(&start,NULL);
        mvdot(a,b,n,n,res1);
        gettimeofday(&end,NULL);
        int noomp = getCPUTime(start,end);
        gettimeofday(&start,NULL);
        mvdot_omp(a,b,n,n,res2);
        gettimeofday(&end,NULL);
        int omp = getCPUTime(start,end);
        for(int i = 0;i < n;i++){
            err[i] = fabs(res1[i] - res2[i]);
        }
        double errnrm = 0.0f;
        dnrm2(err,n,&errnrm);
        printf("n = %d\nError L2nrm = %lf\n",n,errnrm);
        fprintf(bmark,"%d,%d,%d\n",n,omp,noomp);
        for(int i = 0;i < n;i++){
            free(a[i]);
        }
        free(b);
        free(res1);
        free(res2);
        free(err);
    }
    return 0;
}