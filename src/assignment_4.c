#include <fdm.h>
#include <iter_solv.h>
#include <stdio.h>
#include <mpi.h>

int main(int argc, char** argv){
    MPI_Init(&argc, &argv);
    int nproc;
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    return 0;
}