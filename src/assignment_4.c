#include <fdm.h>
#include <iter_solv.h>
#include <stdio.h>
#include <mpi.h>

int main(int argc, char** argv){
    MPI_Init(&argc, &argv);
    int nproc;
    int rank;
    int Ng, N, M;
    int iter = 0;
    double delta;
    double *up, *down, *left, *right;
    double *A, *b, *x;
    double upg, downg, leftg, rightg;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // Initialization
    if(rank == 0){
        printf("Enter the number of points: \n");
        scanf("%d",&Ng);
        printf("Enter the global boundary conditions: \n");
        printf("1) u(0,y) [for all 0 < y < 1] = ");
        scanf("%lf",&leftg);
        printf("2) u(x,1) [for all 0 <= x <= 1] = ");
        scanf("%lf",&upg);
        printf("3) u(1,y) [for all 0 < y < 1] = ");
        scanf("%lf",&rightg);
        printf("4) u(x,0) [for all 0 <= x <= 1] = ");
        scanf("%lf",&downg);
        M = (int)sqrt((double)Ng);
        delta = 1.0/(M - 1);
        MPI_Bcast(&Ng, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&delta, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&M, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Bcast(&upg, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&downg, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&leftg, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&rightg, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
    if(rank < Ng%nproc){
        N = Ng/nproc + 1;
    }
    else{
        N = Ng/nproc;
    }
    if(rank == 0 || rank == nproc - 1)
        N = N + 1; // Leftmost overlaps only on the right side and rightmost overlaps only on the left side
    else
        N = N + 2; // Left and right sides overlap, hence must add two new columns (overlapped points)
    up = (double*) malloc((N-2) * sizeof(double));
    down = (double*) malloc((N-2) * sizeof(double));
    left = (double*) malloc((M-2) * sizeof(double));
    right = (double*) malloc((M-2) * sizeof(double));
    for(int i = 0;i < N;i++){
        up[i] = upg;
        down[i] = downg;
    }
    if(rank == 0){
        for(int i = 0;i < M;i++){
            left[i] = leftg;
        }
    }
    if(rank == nproc - 1){
        for(int i = 0;i < M;i++){
            right[i] = rightg;
        }
    }
    
    // Setting boundary conditions and init
    if(iter == 0){
        // Initial guess
        x = (double*) malloc((N-2)*(M-2) * sizeof(double));
        for(int i = 0;i < (N-2)*(M-2);i++){
            x[i] = 1.0;
        }
        if(rank == 0){
            for(int i = 0;i < N-2;i++){
                right[i] = 1.0;
            }
        }
        if(rank == nproc - 1){
            for(int i = 0;i < N-2;i++){
                left[i] = 1.0;
            }
        }
        if(rank > 0 && rank < nproc - 1){
            for(int i = 0;i < N-2;i++){
                left[i] = 1.0;
                right[i] = 1.0;
            }
        }
    }
    else{
        // Update boundary conditions from neighboring processes
        if(rank == 0){
            
        }
        if(rank == nproc - 1){
            
        }
        if(rank > 0 && rank < nproc - 1){
            
        }
    }

    // Solver
    generateMat(N + 1, M, down, up, left, right, delta, delta, A, b);

    MPI_Finalize();
    return 0;
}