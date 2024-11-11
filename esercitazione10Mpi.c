#include <stdio.h>
#include <mpi.h>
#include <sys/time.h>

double get_cur_time();
void laplace(float *, float *, float *, float *, int, int, int);

int main(int argc, char **argv)
{
    double t1, t2;

    int nproc, myid, prev, next;
    int N, i, j, ifirst, iter, Niter, LD;
    float *A, *Anew, *daprev, *danext;
    MPI_Status status;

    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    printf("hello from %d di %d processi \n", myid, nproc);
    sleep(1);

    N = 400;
    Niter = 8000;
    LD = 500;
    A = (float *)malloc(500 * 500 * sizeof(float));
    Anew = (float *)malloc(500 * 500 * sizeof(float));
    daprev = (float *)malloc(500 * sizeof(float));
    danext = (float *)malloc(500 * sizeof(float));

    // inizializzazione matrice

    for (i = 0; i < N / nproc; i++)
    { // tutta la matrice locale = 0
        for (j = 0; j < N; j++)
        {
            A[i * LD + j] = 0.;
        }
    }
    if (myid == 0)
        for (j = 0; j < N; j++)
            A[0 * LD + j] = j; // prima riga matrice del proc id=0  da 0 a 390

    if (myid == nproc - 1)
        for (j = 0; j < N; j++)
            A[(N / nproc - 1) * LD + j] = N - 1 - j; // ultima riga matrice del proc id=nproc-1 da 390 a 0

    ifirst = myid * N / nproc;

    for (i = 0; i < N / nproc; i++)
    {
        A[i * LD + 0] = ifirst + i;                // bordo sinistro da ifirst a ilast-1 in ogni proc
        A[i * LD + N - 1] = N - 1 - A[i * LD + 0]; // A[i][0] + A[i][N-1] = 0 sempre
    }

    if (myid == 0)
        printf("\n esecuzione con N = %d  e %d iterazioni\n\n", N, Niter);

    t1 = get_cur_time();

    laplace(A, Anew, daprev, danext, N, LD, Niter);

    t2 = get_cur_time();

    if (myid == 0)
        printf("con %d processi, il tempo e' %f\n", nproc, t2 - t1);

    sleep(1);
    if (myid == 0)
        printf("prima  %d -->   %f  %f  \n", myid, A[1 * LD + 1], A[1 * LD + 398]);
    if (myid == 3)
        printf("centro %d -->   %f  %f  \n", myid, A[49 * LD + 199], A[49 * LD + 200]);
    if (myid == 4)
        printf("centro %d -->   %f  %f  \n", myid, A[00 * LD + 199], A[00 * LD + 200]);
    if (myid == 7)
        printf("ultima %d -->   %f  %f  \n", myid, A[48 * LD + 1], A[48 * LD + 398]);

    MPI_Finalize();
}

void laplace(float *A, float *Anew, float *daprev, float *danext, int N, int LD, int Niter)
{
    int iter, i, j;
    int myid, nproc;
    MPI_Status status;

    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    for (iter = 0; iter < Niter; iter++)
    {
        // Se id != 0, invia la prima riga al processo precedente e ricevi l'ultima riga da questo
        if (myid != 0)
        {
            MPI_Send(A, LD, MPI_FLOAT, myid - 1, 0, MPI_COMM_WORLD);
            MPI_Recv(daprev, LD, MPI_FLOAT, myid - 1, 0, MPI_COMM_WORLD, &status);
        }

        // Se id != NP-1, invia l'ultima riga al processo successivo e ricevi la prima riga da questo
        if (myid != nproc - 1)
        {
            MPI_Send(A + (N / nproc - 1) * LD, LD, MPI_FLOAT, myid + 1, 0, MPI_COMM_WORLD);
            MPI_Recv(danext, LD, MPI_FLOAT, myid + 1, 0, MPI_COMM_WORLD, &status);
        }

        // Calcola la parte interna della matrice B
        for (i = 1; i < N / nproc - 1; i++)
        {
            for (j = 1; j < N - 1; j++)
            {
                Anew[i * LD + j] = 0.25 * (A[(i - 1) * LD + j] + A[(i + 1) * LD + j] +
                                           A[i * LD + (j - 1)] + A[i * LD + (j + 1)]);
            }
        }

        // Calcola la prima riga di B se id != 0
        if (myid != 0)
        {
            for (j = 1; j < N - 1; j++)
            {
                Anew[0 * LD + j] = 0.25 * (daprev[j] + A[1 * LD + j] +
                                           A[0 * LD + (j - 1)] + A[0 * LD + (j + 1)]);
            }
        }

        // Calcola l'ultima riga di B se id != NP-1
        if (myid != nproc - 1)
        {
            for (j = 1; j < N - 1; j++)
            {
                Anew[(N / nproc - 1) * LD + j] = 0.25 * (A[(N / nproc - 2) * LD + j] + danext[j] +
                                                         A[(N / nproc - 1) * LD + (j - 1)] + A[(N / nproc - 1) * LD + (j + 1)]);
            }
        }

        // Copia la parte interna di B in A
        for (i = 1; i < N / nproc - 1; i++)
        {
            for (j = 1; j < N - 1; j++)
            {
                A[i * LD + j] = Anew[i * LD + j];
            }
        }

        // Copia la prima riga di B in A se id != 0
        if (myid != 0)
        {
            for (j = 1; j < N - 1; j++)
            {
                A[0 * LD + j] = Anew[0 * LD + j];
            }
        }

        // Copia l'ultima riga di B in A se id != NP-1
        if (myid != nproc - 1)
        {
            for (j = 1; j < N - 1; j++)
            {
                A[(N / nproc - 1) * LD + j] = Anew[(N / nproc - 1) * LD + j];
            }
        }
    }
}

void laplace_nb(float *A, float *Anew, float *daprev, float *danext, int N, int LD, int Niter)
{
    int iter, i, j;
    int myid, nproc;
    MPI_Status status;
    MPI_Request request1, request2, request3, request4;

    

    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);

    for (iter = 0; iter < Niter; iter++)
    {
        // Se id != 0, invia la prima riga al processo precedente e ricevi l'ultima riga da questo
        if (myid != 0)
        {
            MPI_Isend(A, LD, MPI_FLOAT, myid - 1, 0, MPI_COMM_WORLD, &request1);
            MPI_Irecv(daprev, LD, MPI_FLOAT, myid - 1, 0, MPI_COMM_WORLD, &request2);
        }

        

        // Se id != NP-1, invia l'ultima riga al processo successivo e ricevi la prima riga da questo
        if (myid != nproc - 1)
        {
            MPI_Isend(A + (N / nproc - 1) * LD, LD, MPI_FLOAT, myid + 1, 0, MPI_COMM_WORLD, &request3);
            MPI_Irecv(danext, LD, MPI_FLOAT, myid + 1, 0, MPI_COMM_WORLD, &request4);
        }

        

        // Calcola la parte interna della matrice B
        for (i = 1; i < N / nproc - 1; i++)
        {
            for (j = 1; j < N - 1; j++)
            {
                Anew[i * LD + j] = 0.25 * (A[(i - 1) * LD + j] + A[(i + 1) * LD + j] +
                                           A[i * LD + (j - 1)] + A[i * LD + (j + 1)]);
            }
        }

        MPI_Wait(&request1, &status);
        MPI_Wait(&request2, &status);


        // Calcola la prima riga di B se id != 0
        if (myid != 0)
        {
            for (j = 1; j < N - 1; j++)
            {
                Anew[0 * LD + j] = 0.25 * (daprev[j] + A[1 * LD + j] +
                                           A[0 * LD + (j - 1)] + A[0 * LD + (j + 1)]);
            }
        }

           MPI_Wait(&request2, &status);

        // Calcola l'ultima riga di B se id != NP-1
        if (myid != nproc - 1)
        {
            for (j = 1; j < N - 1; j++)
            {
                Anew[(N / nproc - 1) * LD + j] = 0.25 * (A[(N / nproc - 2) * LD + j] + danext[j] +
                                                         A[(N / nproc - 1) * LD + (j - 1)] + A[(N / nproc - 1) * LD + (j + 1)]);
            }
        }

        // Copia la parte interna di B in A
        for (i = 1; i < N / nproc - 1; i++)
        {
            for (j = 1; j < N - 1; j++)
            {
                A[i * LD + j] = Anew[i * LD + j];
            }
        }

        // Copia la prima riga di B in A se id != 0
        if (myid != 0)
        {
            for (j = 1; j < N - 1; j++)
            {
                A[0 * LD + j] = Anew[0 * LD + j];
            }
        }

        // Copia l'ultima riga di B in A se id != NP-1
        if (myid != nproc - 1)
        {
            for (j = 1; j < N - 1; j++)
            {
                A[(N / nproc - 1) * LD + j] = Anew[(N / nproc - 1) * LD + j];
            }
        }
    }
}
double get_cur_time()
{
    struct timeval tv;
    // struct timezone  tz;
    double cur_time;

    gettimeofday(&tv, NULL);
    cur_time = tv.tv_sec + tv.tv_usec / 1000000.0;

    return cur_time;
}
