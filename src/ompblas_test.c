#include <ompblas.h>
#include <cblas.h>
#include "sys/time.h"
#include "time.h"
#include "unistd.h"
double randd(double min, double max){
    return min + ((double)rand()/(double)RAND_MAX)*(max-min);
}
int main(){
    int n;
    printf("Enter n = ");
    scanf("%d",&n);
    printf("\n");
    double min,max;
    min = -1.0;
    max = 1.0;
    double *a, *b, res;
    a = (double*)malloc(n*sizeof(double)); 
    b = (double*)malloc(n*sizeof(double)); 
    for(int i = 0;i < n;i++){
        a[i] = randd(min,max);
        b[i] = randd(min,max);
    }
    set_num_threads(2);
    // struct timeval start, end;
    // gettimeofday(&start,NULL)
    printVect(a,n);
    printVect(b,n);
    vvdot(a,b,n,&res);

    printf("Result (Using ompblas) = %lf\n",res);
    res = cblas_ddot(n,a,1,b,1);
    printf("Result (Using cblas) = %lf\n",res);
    return 0;
}