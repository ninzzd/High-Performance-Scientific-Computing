#include "iter_solv.h"
#include "opt_solv.h"
#include "fdm.h"
/*#include "sys/time.h"
#include "time.h"
#include "unistd.h"
*/
int main(){
    double **a, *b, *x0, *x;
    int n;
    FILE *kmat, *fvec, *kinfo;
    fGetMat(&kmat,&fvec,&kinfo);
    getMat(&a,&b,&n,kmat,fvec,kinfo);
    fclose(kmat);
    fclose(fvec);
    fclose(kinfo);
    
    x0 = (double*)malloc(n*sizeof(double));
    x = (double*)malloc(n*sizeof(double));
    for(int i = 0;i < n;i++){
        x0[i] = 0.0f;
    }
    int iter_bicgstab = bicgstab(a,b,x0,n,0.000001,0,x);
    printf("Number of iterations: %d\n",iter_bicgstab);
    printf("Solution:\n");
    printVect(x,n);
    free(b);
    free(x0);
    free(x);
    for(int i = 0;i < n;i++){
        free(a[i]);
    }
    return 0;
}