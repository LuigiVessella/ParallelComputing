#include <stdio.h>
#include <mpi.h>

int main(int argc, char *argv[])
{

    int nproc, myid;

    MPI_Status status;

    /* tipo MPI definito in mpi.h */

    float a[2];

    /* array da inviare
     */

    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    if (myid == 1)
    {
        a[0] = 1, a[1] = 2;
        MPI_Send(a, 2, MPI_FLOAT,
                 0, 10, MPI_COMM_WORLD); // se sei il processo con id 1 manda al processo 0 l'array. Il processo 0 lo stamperà
    }
    else
    {
        MPI_Recv(a, 2, MPI_FLOAT, 1, 10, MPI_COMM_WORLD, &status);
        printf("%d: a[0]=%f a[1]=%f\n", myid, a[0], a[1]);
    }
    MPI_Finalize();
}
