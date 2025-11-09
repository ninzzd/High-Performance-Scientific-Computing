#include <iter_solv.h>
void printMat(double** a, double* b, double* x, int n){
    printf("\nMatrix A:\n");
    for(int i = 0;i < n;i++){
        for(int j = 0;j < n;j++){
            printf("%lf ",a[i][j]);
        }
        printf("\n");
    }
    printf("\nMatrix B:\n");
    for(int i = 0;i < n;i++){
        printf("%lf\n",b[i]);
    }
    printf("\nMatrix X:\n");
    for(int i = 0;i < n;i++){
        printf("%lf\n",x[i]);
    }
}
// void printVect(double* x, int n){
//     printf("\n[");
//     for(int i = 0;i < n;i++){
//         printf("%lf ",x[i]);
//     }
//     printf("]\n");
// }
void printIter(long long int iter){
    const char* suffix;
        switch(iter%10){
            case 1:
                suffix = "st";
                break;
            case 2:
                suffix = "nd";
                break;
            case 3:
                suffix = "rd";
                break;
            default:
                suffix = "th";
        }
    printf("\n%d%s iteration:\n",iter,suffix);
}
// mode = 0 : Non-verbose, mode = 1 : Verbose
long long int jacobiSolver(double** a, double* b, double* x_0, int n, double eta, int mode, double* result){
    double* x_k = (double*)malloc(n*sizeof(double));
    double* x_k1 = (double*)malloc(n*sizeof(double));
    long long int iter = 0;
    double max_err;
    do{
        if(iter == 0)
            memcpy(x_k,x_0,(size_t)n*sizeof(double));
        else
            memcpy(x_k,x_k1,(size_t)n*sizeof(double));
        for(int i = 0;i < n;i++){
            double sum = b[i];
            for(int j = 0;j < n; j++){
                if(j != i)
                    sum -= a[i][j]*x_k[j];
            }
            x_k1[i] = sum/a[i][i];
        }
        iter++;
        max_err = 0.0;
        for(int i = 0;i < n;i++){
            if(fabs(x_k1[i] - x_k[i]) > max_err)
                max_err = fabs(x_k1[i] - x_k[i]);
        }
        if(mode){
            printIter(iter);
            printVect(x_k1,n);
            printf("Max Error = %lf\n",max_err);
        }
    }while(max_err > eta);
    memcpy(result,x_k1,(size_t)n*sizeof(double));
    free(x_k);
    free(x_k1);
    return iter;
}
long long int gaussSeidelSolver(double** a, double* b, double* x_0, int n, double eta, int mode, double* result){
    double* x_k = (double*)malloc(n*sizeof(double));
    double* x_k1 = (double*)malloc(n*sizeof(double));
    long long int iter = 0;
    double max_err;
    do{
        if(iter == 0)
            memcpy(x_k,x_0,(size_t)n*sizeof(double));
        else
            memcpy(x_k,x_k1,(size_t)n*sizeof(double));
        for(int i = 0;i < n;i++){
            double sum = b[i];
            for(int j = 0;j < n; j++){
                if(j < i)
                    sum -= a[i][j]*x_k1[j];
                else if(j > i)
                    sum -= a[i][j]*x_k[j];
            }
            x_k1[i] = sum/a[i][i];
        }
        iter++;
        max_err = 0.0;
        for(int i = 0;i < n;i++){
            if(fabs(x_k1[i] - x_k[i]) > max_err)
                max_err = fabs(x_k1[i] - x_k[i]);
        }
        if(mode){
            printIter(iter);
            printVect(x_k1,n);
            printf("Max Error = %lf\n",max_err);
        }
    }while(max_err > eta);
    memcpy(result,x_k1,(size_t)n*sizeof(double));
    free(x_k);
    free(x_k1);
    return iter;
}
long long int sorSolver(double** a, double* b, double* x_0, int n, double eta, double w, int mode, double* result){
    double* x_k = (double*)malloc(n*sizeof(double));
    double* x_k1 = (double*)malloc(n*sizeof(double));
    long long int iter = 0;
    double max_err;
    do{
        if(iter == 0)
            memcpy(x_k,x_0,(size_t)n*sizeof(double));
        else
            memcpy(x_k,x_k1,(size_t)n*sizeof(double));
        for(int i = 0;i < n;i++){
            double sum = b[i];
            for(int j = 0;j < n; j++){
                if(j < i)
                    sum -= a[i][j]*x_k1[j];
                else if(j > i)
                    sum -= a[i][j]*x_k[j];
            }
            x_k1[i] = sum/a[i][i];
            x_k1[i] = x_k[i] + w*(x_k1[i] - x_k[i]);
        }
        iter++;
        double* err = (double*)malloc(n*sizeof(double));
        for(int i = 0;i < n;i++)
            err[i] = x_k1[i] - x_k[i];
        max_err = cblas_dnrm2(n,err,1);
        // for(int i = 0;i < n;i++){
        //     if(fabsf(x_k1[i] - x_k[i]) > max_err)
        //         max_err = fabsf(x_k1[i] - x_k[i]);
        // }

        if(mode){
            printIter(iter);
            printVect(x_k1,n);
            printf("\n Error vector:\n");
            for(int i = 0;i < n;i++){
                printf("%lf ",x_k1[i] - x_k[i]);
            }
            printf("\n");
            printf("Max Error = %lf\n",max_err);
        }
    }while(max_err > eta);
    memcpy(result,x_k1,(size_t)n*sizeof(double));
    free(x_k);
    free(x_k1);
    return iter;
}
double sorOpt(double **a, int n){
    double* lambda_r = (double*)malloc(n*sizeof(double));
    double* lambda_i = (double*)malloc(n*sizeof(double));
    double* a_flat = (double*)malloc(n*n*sizeof(double));
    double diag_prod = 1.0;
    for(int i = 0;i < n;i++){
        for(int j = 0;j < n;j++){
            a_flat[i*n + j] = a[i][j];
            // if(i == j)
            //     diag_prod *= a[i][j];
        }
    }
    double* d_flat = (double*)malloc(n*n*sizeof(double));
    for(int i = 0;i < n;i++){
        for(int j = 0;j < n;j++){
            if(i == j)
                d_flat[i*n + j] = 1.0/a[i][j];
            else
                d_flat[i*n + j] = 0.0;
        }
    }
    double* e_f_flat = (double*)malloc(n*n*sizeof(double));
    for(int i = 0;i < n;i++){
        for(int j = 0;j < n;j++){
            if(i == j)
                e_f_flat[i*n + j] = 0.0;
            else
                e_f_flat[i*n + j] = -a[i][j];
        }
    }
    double* g = (double*)malloc(n*n*sizeof(double));
    memset(g,0,n*n*sizeof(double));
    // C <- alpha*A*B + beta*C
    cblas_dgemm(CblasRowMajor,CblasNoTrans,CblasNoTrans,n,n,n,1.0,d_flat,n,e_f_flat,n,0.0,g,n);
    int info = LAPACKE_dgeev(LAPACK_ROW_MAJOR,'N','N',n,g,n,lambda_r,lambda_i,NULL,n,NULL,n);
    double w_opt;
    if(info == 0){
        double sr = -__DBL_MAX__;
        for(int i = 0;i < n;i++){
            double mag = sqrt(lambda_r[i]*lambda_r[i] + lambda_i[i]*lambda_i[i]);
            if(mag > sr)
                sr = mag;
        }
        printf("Spectral Radius = %f\n",sr);
        w_opt = 2.0/(1.0 + sqrt(1.0 - sr*sr));
    }
    else{
        printf("Error in computing eigenvalues. Info = %d\n",info);
    }
    free(lambda_r);
    free(lambda_i);
    free(a_flat);
    return w_opt;
}
void fGetMat(FILE** kmat, FILE** fvec, FILE** kinfo){
    char dir[255];
    char kmat_path[255];
    char kinfo_path[255];
    char fvec_path[255];
    printf("Enter the directory (absolute path) containing the files: 1) kinfo.txt 2) kmat.txt 3) Fvec.txt:\n");
    scanf("%s",&dir);
    strcpy(kmat_path,dir);
    strcpy(kinfo_path,dir);
    strcpy(fvec_path,dir);
    strcat(kmat_path,"kmat.txt");
    strcat(kinfo_path,"kinfo.txt");
    strcat(fvec_path,"Fvec.txt");
    *kmat = fopen(kmat_path, "r");
    *kinfo = fopen(kinfo_path, "r");
    *fvec = fopen(fvec_path, "r");
}
void getMat(double*** a, double** b, int* n, FILE* kmat, FILE* fvec, FILE* kinfo){
    fscanf(kinfo,"%d",n);
    *a = (double**)malloc(*n*sizeof(double*));
    for(int i = 0;i < *n;i++){
        (*a)[i] = (double*)malloc(*n*sizeof(double));
    }
    *b = (double*)malloc(*n*sizeof(double));
    for(int i = 0;i < *n;i++){
        for(int j = 0;j < *n;j++){
            fscanf(kmat,"%lf",&(*a)[i][j]);
        }
    }
    for(int i = 0;i < *n;i++){
        fscanf(fvec,"%lf",&(*b)[i]);
    }
}