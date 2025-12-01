#include <fdm.h>
#include <iter_solv.h>
#include <stdio.h>
#include <mpi/mpi.h>

int main(int argc, char** argv){
    MPI_Init(&argc, &argv);
    int nproc;
    int rank;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int Ng, N, M;
    int iter, reqLen;
    double delta;
    double *up, *down, *left, *right; 
    double *sendRight, *sendLeft;
    double *A, *b, *x, *x_next;
    double upg, downg, leftg, rightg;
    double err, errg, errg_threshold = 1e-6;
    MPI_Request *reqs;
    iter = 1;
    err = -__DBL_MAX__;
    // Initialization
    if(rank == 0){
        printf("Enter the number of points: \n");
        scanf("%d",&Ng);
        printf("Enter the global boundary conditions: \n");
        printf("1) u(0,y) [for all 0 < y < 1] [LEFT] = \n");
        scanf("%lf",&leftg);
        printf("2) u(x,1) [for all 0 <= x <= 1] [UP] = \n");
        scanf("%lf",&upg);
        printf("3) u(1,y) [for all 0 < y < 1] [RIGHT] = \n");
        scanf("%lf",&rightg);
        printf("4) u(x,0) [for all 0 <= x <= 1] [DOWN] = \n");
        scanf("%lf",&downg);
        M = (int)sqrt((double)Ng);
        delta = 1.0/(M - 1);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&Ng, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&delta, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&M, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&upg, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&downg, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&leftg, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&rightg, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    printf("Rank = %d, Process ID = %d\n", rank, getpid());
    // Load balancing (Expected: Ng >> nproc)
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
    
    A = (double*)malloc((N-2) * (M-2) * (N-2) * (M-2) * sizeof(double));
    b = (double*)malloc((N-2) * (M-2) * sizeof(double));
    up = (double*) malloc((N-2) * sizeof(double));
    down = (double*) malloc((N-2) * sizeof(double));
    left = (double*) malloc((M-2) * sizeof(double));
    right = (double*) malloc((M-2) * sizeof(double));
    x = (double*) malloc((N-2)*(M-2) * sizeof(double));
    x_next = (double*) malloc((N-2)*(M-2) * sizeof(double));
    sendRight = malloc((M-2)*sizeof(double));
    sendLeft = malloc((M-2)*sizeof(double));
    reqLen = (rank == 0 || rank == nproc - 1) ? 2 : 4;
    reqs = (MPI_Request*) malloc(reqLen*sizeof(MPI_Request));
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
    do{
        // Setting boundary conditions and init
        if(iter == 1){
            // Initial guess
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
                    right[i] = 1.0;
                    left[i] = 1.0;
                }
            }
        }
        else{
            // Update boundary conditions from neighboring processes
            int reqCount = 0;
            // Left to right communication
            if(rank < nproc - 1){ // Initializing left to right send
                for(int i = 0;i < M-2; i++){
                    sendRight[i] = x[(N-3) + i*(N-2)];
                }
                MPI_Isend(sendRight, M-2, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &reqs[reqCount++]); // Left to right send
                MPI_Irecv(right, M-2, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD, &reqs[reqCount++]); // Right to left recv
            }
            if(rank > 0){
                for(int i = 0;i < M-2; i++){
                    sendLeft[i] = x[i*(N-2)];
                }
                MPI_Isend(sendLeft, M-2, MPI_DOUBLE, rank - 1, 1, MPI_COMM_WORLD, &reqs[reqCount++]); // Right to left send
                MPI_Irecv(left, M-2, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &reqs[reqCount++]); // Left to right recv
            }
        }
        MPI_Waitall(reqLen, reqs, MPI_STATUSES_IGNORE);
        // Solver
        generateMat(N, M, down, up, left, right, delta, delta, A, b); 
        for(int i = 0;i < (N-2)*(M-2);i++){
            double sum = 0.0;
            for(int j = 0;j < (N-2)*(M-2);j++){
                if(j == i) 
                    continue;
                sum += A[i*((N-2)*(M-2)) + j] * x[j];
            }
            x_next[i] = (b[i] - sum)/A[i*((N-2)*(M-2)) + i];
            err = fmax(err, fabs(x_next[i] - x[i]));
        }
        MPI_Allreduce(&err, &errg, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
        if(rank == 0){
            printf("Global Error after iteration %d: %lf\n", iter, errg);
        }
        MPI_Bcast(&errg, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        memcpy(x, x_next, (N-2)*(M-2)*sizeof(double));
        iter++;
    }while(errg > errg_threshold);
    // Cleanup
    free(reqs);
    free(sendLeft);
    free(sendRight);
    free(up);
    free(down);
    free(left);
    free(right);
    free(A);
    free(b);
    free(x);
    free(x_next);
    MPI_Finalize();
    return 0;
}