#include <stdio.h>
#include <mpi.h>

int main (int argc, char** argv) {

    int nproc, myid, i, A[20000], somma, N;

    MPI_Init(&argc, &argv);
    
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    MPI_Status status;

    N = 16000;
    for(i = 0; i<N; i++){
        A[i] = myid;
    }
    somma = 0;
    for(i = 0; i < nproc; i++) {
        MPI_Send(A, N, MPI_INT, (myid+nproc+1)%nproc, 10, MPI_COMM_WORLD);
        MPI_Recv(A, N, MPI_INT, (myid+nproc-1)%nproc, 10, MPI_COMM_WORLD, &status);
        somma = somma + 1;
    }

    printf("il numero di shift del proc %d e' %d \n", myid, somma);

    MPI_Finalize();
    return 0;
}

