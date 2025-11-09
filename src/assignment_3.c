#include "iter_solv.h"
#include "opt_solv.h"
#include <ompblas.h>
#include "fdm.h"
#include "omp.h"
#include "sys/time.h"
#include "time.h"
#include "unistd.h"

int main(){
    struct timeval start, end;
    double **a, *b, *x0, *x_noomp, *x_omp;
    int n,num_iter = 1000;
    FILE *kmat, *fvec, *kinfo;
    fGetMat(&kmat,&fvec,&kinfo);
    getMat(&a,&b,&n,kmat,fvec,kinfo);
    fclose(kmat);
    fclose(fvec);
    fclose(kinfo);
    
    x0 = (double*)malloc(n*sizeof(double));
    x_noomp = (double*)malloc(n*sizeof(double));
    x_omp = (double*)malloc(n*sizeof(double));
    //err = (double*)malloc(n*sizeof(double));
    for(int i = 0;i < n;i++){
        x0[i] = 0.0f;
    }
    double noomp_avg,omp_avg;
    omp_avg = noomp_avg = 0.0;
    for(int j = 2;j <= 8;j*=2){
        set_num_threads(j);
        for(int iter = 1;iter <= num_iter;iter++){
            gettimeofday(&start,NULL);
            bicgstab(a,b,x0,n,0.000001,0,x_noomp);
            gettimeofday(&end,NULL);
            int t_noomp = getCPUTime(start,end,0);
            noomp_avg += (double)t_noomp;

            gettimeofday(&start,NULL);
            bicgstabOMP(a,b,x0,n,0.000001,0,x_omp);
            gettimeofday(&end,NULL);
            int t_omp = getCPUTime(start,end,0);
            omp_avg += (double)t_omp;
        }
        noomp_avg /= (double)num_iter;
        omp_avg /= (double)num_iter;
        // for(int i = 0;i < n;i++){
        //     err[i] = fabs(x_omp[i] - x_noomp[i]);
        // }
        // double errnrm2;
        // dnrm2(err,n,&errnrm2);
        printf("Number of Nodes = %d\n",n);
        printf("Number of Threads = %d\n",get_num_threads());
        // printf("Number of iterations (No-OMP) = %d\n",iter_bicgstab_noomp);
        // printf("Time taken (No-OMP) = %d\n",t_noomp);
        // printf("Number of iterations (OMP): %d\n",iter_bicgstab_omp);
        // printf("Time taken (OMP) = %d\n",t_omp);
        // printf("Error L2 NRM = %lf\n",errnrm2);
        printf("Avg Time (No-OMP) (%d Iterations) = %lf\n",num_iter,noomp_avg);
        printf("Avg Time (OMP) (%d Iterations) = %lf\n",num_iter,omp_avg);
    }
    free(b);
    free(x0);
    free(x_noomp);
    free(x_omp);
    for(int i = 0;i < n;i++){
        free(a[i]);
    }
    // printf("Maximum thread count limit: %d\n",omp_get_max_threads());
    return 0;
}