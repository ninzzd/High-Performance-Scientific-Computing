#include "opt_solv.h"
int steepestDescent(double **a, double *b, double* x0, int n, double eta, int mode, double* result){
    double *xk = (double*)malloc(n*sizeof(double));
    double *xk1 = (double*)malloc(n*sizeof(double));
    double *rk = (double*)malloc(n*sizeof(double));
    double rkl2nrm;
    int iter = 0;
    do{
        iter++;
        if(iter == 1)
            memcpy(xk,x0,n*sizeof(double));
        else
            memcpy(xk,xk1,n*sizeof(double));
        for(int i = 0;i < n;i++){
            double ax = 0.0;
            for(int j = 0;j < n;j++){
                ax += a[i][j]*xk[j];
            }
            rk[i] = b[i] - ax;
        }
        double alpha = 0.0;
        double rl2sq = 0.0;
        double rtar = 0.0;
        for(int i = 0;i < n;i++){
            rl2sq += rk[i]*rk[i];
        }
        for(int i = 0;i < n;i++){
            for(int j = 0;j < n;j++){
                rtar += rk[i]*a[i][j]*rk[j];
            }
        }
        alpha = rl2sq/rtar;
        for(int i = 0;i < n;i++){
            xk1[i] = xk[i] + alpha*rk[i];
        }   
        rkl2nrm = cblas_dnrm2(n,rk,1);
        if(mode){
            printIter(iter);
            printf("x(%d) = ",iter);
            printVect(xk1,n);
            printf("L2 Norm of r(k) = %lf\n",rkl2nrm);
        }
        iter++;
    }while(rkl2nrm > eta);
    memcpy(result,xk1,n*sizeof(double));
    free(xk);
    free(xk1);
    free(rk);
    return iter;
}
int minimalResidual(double **a, double *b, double* x0, int n, double eta, int mode, double* result){
    double *xk = (double*)malloc(n*sizeof(double));
    double *rk = (double*)malloc(n*sizeof(double));
    double *pk = (double*)malloc(n*sizeof(double));
    double err;
    int iter = 0;
    for(int i = 0;i < n;i++){
        double ax = 0.0;
        for(int j = 0;j < n;j++){
            ax += a[i][j]*x0[j];
        }
        rk[i] = b[i] - ax;
    }
    do{
        iter++;
        if(iter == 1){
            memcpy(xk,x0,n*sizeof(double));
        }
        for(int i = 0;i < n;i++){
            double row = 0.0;
            for(int j = 0;j < n;j++){
                row += a[i][j]*rk[j];
            }
            pk[i] = row;
        }
        double alpha = cblas_ddot(n,rk,1,pk,1)/pow(cblas_dnrm2(n,pk,1),2.0);
        for(int i = 0;i < n;i++){
            xk[i] += alpha*rk[i];
        }
        for(int i = 0;i < n;i++){
            rk[i] -= alpha*pk[i];
        }
        err = cblas_dnrm2(n,rk,1);
    }while(err > eta);
    memcpy(result,xk,n*sizeof(double));
    return iter;
    free(xk);
    free(rk);
    free(pk);
    return 0;
}
int conjugateGradient(double **a, double *b, double *x0, int n, double eta, int mode, double *result){
    double *xk = (double*)malloc(n*sizeof(double));
    double *rk = (double*)malloc(n*sizeof(double));
    double *pk = (double*)malloc(n*sizeof(double));
    double *apk = (double*)malloc(n*sizeof(double));
    double err;
    int iter = 0;
    for(int i = 0;i < n;i++){
        double ax = 0.0;
        for(int j = 0;j < n;j++){
            ax += a[i][j]*x0[j];
        }
        rk[i] = b[i] - ax;
        pk[i] = rk[i];
    }
    do{
        iter++;
        if(iter == 1){
            memcpy(xk,x0,n*sizeof(double));
        }
        double *apk = (double*)malloc(n*sizeof(double));
        for(int i = 0;i < n;i++){
            double row = 0.0;
            for(int j = 0;j < n;j++){
                row += a[i][j]*pk[j];
            }
            apk[i] = row;
        }
        double alpha = cblas_ddot(n,rk,1,rk,1)/cblas_ddot(n,pk,1,apk,1);
        for(int i = 0;i < n;i++){
            xk[i] += alpha*pk[i];
        }
        double rt_dot_r = cblas_ddot(n,rk,1,rk,1);
        for(int i = 0;i < n;i++){
            rk[i] -= alpha*apk[i];
        }
        double beta = cblas_ddot(n,rk,1,rk,1)/rt_dot_r;
        for(int i = 0;i < n;i++){
            pk[i] = rk[i] + beta*pk[i];
        }
        err = cblas_dnrm2(n,rk,1);
        free(apk);
    }while(err > eta);
    memcpy(result,xk,n*sizeof(double));
    free(xk);
    free(rk);
    free(pk);
    free(apk);
    return iter;
}
int bicgstab(double **a, double *b, double *x0, int n, double eta, int mode, double *result){
    double *xk = (double*)malloc(n*sizeof(double));
    double *rk = (double*)malloc(n*sizeof(double));
    double *r0_ = (double*)malloc(n*sizeof(double));
    double *pk = (double*)malloc(n*sizeof(double));
    double err;
    memcpy(xk,x0,n*sizeof(double));
    // Initializing rk to r0, pk = p0 = r0, r0* = random (arbitrary)
    for(int i = 0;i < n;i++){
        double ax = 0.0;
        for(int j = 0;j < n;j++){
            ax += a[i][j]*x0[j];
        }
        rk[i] = b[i] - ax;
        pk[i] = rk[i];
        r0_[i] = (double)rand()/(double)RAND_MAX; // Is rand() costly?
    }
    int iter = 0;
    do{
        iter++;
        double *ap = (double*)malloc(n*sizeof(double));
        double *s = (double*)malloc(n*sizeof(double));
        mvdot(a,pk,n,n,ap);
        double rk_r0,ap_r0;
        vvdot(rk,r0_,n,&rk_r0);
        vvdot(ap,r0_,n,&ap_r0);
        double alpha = rk_r0/ap_r0;
        for(int i = 0;i < n;i++){
            s[i] = rk[i] - alpha*ap[i];
        }
        double *as = (double*)malloc(n*sizeof(double));
        mvdot(a,s,n,n,as);
        double as_s,as_as;
        vvdot(as,s,n,&as_s);
        vvdot(as,as,n,&as_as);
        double w = as_s/as_as;
        // double temp = cblas_ddot(n,rk,1,r0_,1); 
        // Replace temp wih rk_r0
        double temp = rk_r0;
        for(int i = 0;i < n;i++){
            xk[i] += alpha*pk[i] + w*s[i];
            rk[i] = s[i] - w*as[i];
        }
        vvdot(rk,r0_,n,&rk_r0);
        double beta = (alpha*rk_r0)/(w*temp);
        for(int i = 0;i < n;i++){
            pk[i] = rk[i] + beta*(pk[i] - w*ap[i]);
        }
        dnrm2(rk,n,&err);
        free(ap);
        free(as);
        free(s);
    }while(err > eta);
    memcpy(result,xk,n*sizeof(double));
    free(xk);
    free(rk);
    free(r0_);
    free(pk);
    return iter;
}
int bicgstabOMP(double **a, double *b, double *x0, int n, double eta, int mode, double *result){
    double *xk = (double*)malloc(n*sizeof(double));
    double *rk = (double*)malloc(n*sizeof(double));
    double *r0_ = (double*)malloc(n*sizeof(double));
    double *pk = (double*)malloc(n*sizeof(double));
    double *ax0 = (double*)malloc(n*sizeof(double));
    double err;
    int i;
    memcpy(xk,x0,n*sizeof(double));
    // Initializing rk to r0, pk = p0 = r0, r0* = random (arbitrary)
    mvdot_omp(a,x0,n,n,ax0);

    #pragma omp parallel for num_threads(n_threads) private(i) default(shared)
    for(i = 0;i < n;i++){
        rk[i] = b[i] - ax0[i];
        pk[i] = rk[i];
        r0_[i] = (double)rand()/(double)RAND_MAX; // Is rand() costly?
    }
    int iter = 0;
    do{
        iter++;
        double *ap = (double*)malloc(n*sizeof(double));
        double *s = (double*)malloc(n*sizeof(double));
        mvdot_omp(a,pk,n,n,ap);
        double rk_r0,ap_r0;
        vvdot(rk,r0_,n,&rk_r0);
        vvdot(ap,r0_,n,&ap_r0);
        double alpha = rk_r0/ap_r0;
        #pragma omp parallel for num_threads(n_threads) firstprivate(alpha) private(i) default(shared)
        for(i = 0;i < n;i++){
            s[i] = rk[i] - alpha*ap[i];
        }

        double *as = (double*)malloc(n*sizeof(double));
        mvdot_omp(a,s,n,n,as);
        double as_s,as_as;
        vvdot(as,s,n,&as_s);
        vvdot(as,as,n,&as_as);
        double w = as_s/as_as;
        // double temp = cblas_ddot(n,rk,1,r0_,1); 
        // Replace temp wih rk_r0
        double temp = rk_r0;

        #pragma omp parallel for num_threads(n_threads) firstprivate(alpha,w) private(i) default(shared)
        for(i = 0;i < n;i++){
            xk[i] += alpha*pk[i] + w*s[i];
            rk[i] = s[i] - w*as[i];
        }

        vvdot(rk,r0_,n,&rk_r0);
        double beta = (alpha*rk_r0)/(w*temp);

        #pragma omp parallel for num_threads(n_threads) firstprivate(beta,w) private(i) default(shared)
        for(i = 0;i < n;i++){
            pk[i] = rk[i] + beta*(pk[i] - w*ap[i]);
        }

        dnrm2(rk,n,&err);
        free(ap);
        free(as);
        free(s);
    }while(err > eta);
    memcpy(result,xk,n*sizeof(double));
    free(xk);
    free(rk);
    free(r0_);
    free(pk);
    return iter;
}
