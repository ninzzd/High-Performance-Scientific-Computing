#include <ompblas.h>
#include "sys/time.h"
#include "time.h"
#include "unistd.h"
double randd(double min, double max){
    return min + ((double)rand()/(double)RAND_MAX)*(max-min);
}
int main(){
    int n;
    double min,max;
    min = -1.0;
    max = 1.0;
    double *a, *b, *res;
    a = (double*)malloc(n*sizeof(double)); 
    b = (double*)malloc(n*sizeof(double)); 
    res = (double*)malloc(n*sizeof(double)); 
    for(int i = 0;i < n;i++){
        a[i] = randd(min,max);
        b[i] = randd(min,max);
    }
    set_num_threads(2);
    // struct timeval start, end;
    // gettimeofday(&start,NULL)
    return 0;
}