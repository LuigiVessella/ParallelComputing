// usare questo programma chiamante per fare test di correttezza per matmatdist
// SOLO CON GLIGLIE DI PROCESSI (NPROW , NPCOL) = (1,1) e (2,2)

#include <unistd.h>
#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>

void matmatikj(int ldA, int ldB, int ldC, double *A, double *B, double *C, int N1, int N2, int N3)
void matmatdist(MPI_Comm, int, int, int, double *, double *, double *, int, int, int, int, int, int, int, int);
void matmatblock(int ldA, int ldB, int ldC, double *A, double *B, double *C, int N1, int N2, int N3, int dbA, int dbB, int dbC);
void matmatthread(int ldA, int ldB, int ldC, double *A, double *B, double *C,int N1, int N2, int N3, int dbA, int dbB, int dbC,int NTrow, int NTcol);
double get_cur_time();

int main(int argc, char *argv[])
{
    int i, j, Nglob, Mglob, Pglob, lda, mcm;
    int dims[2], period[2], coord[2], TROW, TCOL, rank, size;
    int X, Y, Q, R;
    double *A, *B, *C, *D;
    double time1, time2, Ndouble;
   
    MPI_Comm GridCom;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    //
    // qua viene definita la griglia di processi
    // ATTENZIONE: il prodotto dims[0]*dims[1] deve essere uguale al
    // numero di processi lanciati da mpirun nel file.pbs
    //
    dims[0] = 1;
    dims[1] = 2;
    period[0] = 1;
    period[1] = 1;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, period, 0, &GridCom);

    //
    // allocazione dello spazio per i test
    //
    lda = 6144;
    A = (double *)malloc(sizeof(double) * lda * lda);
    B = (double *)malloc(sizeof(double) * lda * lda);
    C = (double *)malloc(sizeof(double) * lda * lda);
    D = (double *)malloc(sizeof(double) * lda * lda);

    // ==================================================
    // test di correttezza risultati. Verificare solo per griglie di processi (1,1) e (2,2)
    // ==================================================

    Nglob = 2;
    Mglob = 4;
    Pglob = 4;
    TROW = 1;
    TCOL = 1;

    MPI_Cart_coords(GridCom, rank, 2, coord);

    // il calcolo del mcm serve solo per la stampa del risultato del test di correttezza
    X = dims[0];
    Y = dims[1];
    while (Y != 0)
    {
        Q = X / Y;
        R = X - Q * Y;
        X = Y;
        Y = R;
    }
    mcm = dims[0] * dims[1] / X;

    //
    // definizione delle matrici di input
    //
    for (i = 0; i < Nglob / dims[0]; i++)
    {
        for (j = 0; j < Mglob / mcm; j++)
        {
            A[i * lda + j] = coord[0] * dims[1] * Mglob / mcm + coord[1] * Mglob / mcm + i * Mglob + j;
        }
    }
    for (i = 0; i < Mglob / mcm; i++)
    {
        for (j = 0; j < Pglob / dims[1]; j++)
        {
            B[i * lda + j] = 10 + coord[0] * dims[1] * Pglob + coord[1] * Pglob / dims[1] + i * Pglob + j;
        }
    }
    for (i = 0; i < Nglob / dims[0]; i++)
    {
        for (j = 0; j < Pglob / dims[1]; j++)
        {
            C[i * lda + j] = 0.0;
        }
    }

    matmatdist(GridCom, lda, lda, lda, A, B, C, Nglob, Mglob, Pglob, 1, 1, 1, TROW, TCOL);

    //
    // stampa delle matrici A, B e C
    //
    for (i = 0; i < Nglob / dims[0]; i++)
    {
        for (j = 0; j < Mglob / mcm; j++)
        {
            printf("MAT A id %d->  %f \n", rank, A[i * lda + j]);
        }
    }
    printf("------------------\n");

    for (i = 0; i < Mglob / mcm; i++)
    {
        for (j = 0; j < Pglob / dims[1]; j++)
        {
            printf("MAT B id %d->  %f \n", rank, B[i * lda + j]);
        }
    }
    printf("------------------\n");
    for (i = 0; i < Nglob / dims[0]; i++)
    {
        for (j = 0; j < Pglob / dims[1]; j++)
        {
            printf("MAT C id %d->  %f \n", rank, C[i * lda + j]);
        }
    }

    // ==================================================
    // test di efficienza
    // ==================================================

    srand(0);
    for (i = 0; i < lda; i++)
    {
        for (j = 0; j < lda; j++)
        {
            *(A + i * lda + j) = (float)rand() / RAND_MAX;
            *(B + i * lda + j) = (float)rand() / RAND_MAX;
            *(C + i * lda + j) = (float)rand() / RAND_MAX;
            *(D + i * lda + j) = *(C + i * lda + j);
        }
    }

    if (rank == 0)
        printf("               N         time       Gflops\n");
    for (Nglob = 2048; Nglob <= 2048 * 3; Nglob = Nglob + 2048)
    {
        Ndouble = Nglob;

        TROW = 1;
        TCOL = 1; // test con 1 thread per processo
        MPI_Barrier(MPI_COMM_WORLD);
        time1 = get_cur_time();
        matmatdist(GridCom, lda, lda, lda, A, B, C, Nglob, Nglob, Nglob, 256, 256, 256, TROW, TCOL);
        time2 = get_cur_time() - time1;
        printf(" proc = %d:   %4d   %4d   %e  %f \n", rank, Nglob, TROW * TCOL, time2, 2 * Ndouble * Ndouble * Ndouble / time2 / 1.e9);

        TROW = 2;
        TCOL = 1; // test con 2 thread per processo
        MPI_Barrier(MPI_COMM_WORLD);
        time1 = get_cur_time();
        matmatdist(GridCom, lda, lda, lda, A, B, C, Nglob, Nglob, Nglob, 256, 256, 256, TROW, TCOL);
        time2 = get_cur_time() - time1;
        printf(" proc = %d:   %4d   %4d   %e  %f \n", rank, Nglob, TROW * TCOL, time2, 2 * Ndouble * Ndouble * Ndouble / time2 / 1.e9);

        TROW = 2;
        TCOL = 2; // test con 4 thread per processo
        MPI_Barrier(MPI_COMM_WORLD);
        time1 = get_cur_time();
        matmatdist(GridCom, lda, lda, lda, A, B, C, Nglob, Nglob, Nglob, 256, 256, 256, TROW, TCOL);
        time2 = get_cur_time() - time1;
        printf(" proc = %d:   %4d   %4d   %e  %f \n", rank, Nglob, TROW * TCOL, time2, 2 * Ndouble * Ndouble * Ndouble / time2 / 1.e9);

        TROW = 4;
        TCOL = 2; // test con 4 thread per processo
        MPI_Barrier(MPI_COMM_WORLD);
        time1 = get_cur_time();
        matmatdist(GridCom, lda, lda, lda, A, B, C, Nglob, Nglob, Nglob, 256, 256, 256, TROW, TCOL);
        time2 = get_cur_time() - time1;
        printf(" proc = %d:   %4d   %4d   %e  %f \n", rank, Nglob, TROW * TCOL, time2, 2 * Ndouble * Ndouble * Ndouble / time2 / 1.e9);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
}

void matmatikj(int ldA, int ldB, int ldC, double *A, double *B, double *C, int N1, int N2, int N3)
{

    int i, j, k;
    for (i = 0; i < N1; i++)
    {
        for (k = 0; k < N2; k++)
        {
            for (j = 0; j < N3; j++)
            {
                C[i * ldC + j] = C[i * ldC + j] + (A[i * ldA + k] * B[k * ldB + j]);
            }
        }
    }
}

void matmatblock(int ldA, int ldB, int ldC, double *A, double *B, double *C, int N1, int N2, int N3, int dbA, int dbB, int dbC)
{
    int i, j, k;
    int ii, jj, kk;

    printf("ciao sono qui\n");

    for (ii = 0; ii < N1; ii += dbA)
    {
        for (kk = 0; kk < N2; kk += dbB)
        {
            for (jj = 0; jj < N3; jj += dbC)
            {
                matmatikj(ldA, ldB, ldC, &A[ii * ldA + kk], &B[kk * ldB + jj], &C[ii * ldC + jj], dbA, dbB, dbC);
            }
        }
    }
}

void matmatthread(int ldA, int ldB, int ldC, double *A, double *B, double *C,
                  int N1, int N2, int N3, int dbA, int dbB, int dbC,
                  int NTrow, int NTcol)
{

    int thread_id, IDi, IDj, start_i, end_i, start_j, end_j, block_rows, block_cols;
    omp_set_num_threads(NTrow * NTcol);
#pragma omp parallel private(thread_id, IDi, IDj, start_i, end_i, start_j, end_j)
    {

        thread_id = omp_get_thread_num();
        IDi = thread_id / NTcol; // Indice della riga
        IDj = thread_id % NTcol; // Indice della colonna

        // Calcolare i range del blocco gestito dal thread
        start_i = IDi * (N1 / NTrow);
        end_i = (IDi + 1) * (N1 / NTrow);

        start_j = IDj * (N3 / NTcol);
        end_j = (IDj + 1) * (N3 / NTcol);

        // Determinare le dimensioni dei blocchi da calcolare
        block_rows = end_i - start_i;
        block_cols = end_j - start_j;

        printf("processo %d%d calcola %d %d della mat\n", IDi, IDj, block_rows, block_cols);

        // Chiamata alla funzione `matmatblock` per calcolare il blocco

        matmatblock(ldA, ldB, ldC,
                    &A[start_i * ldA],           // Offset riga del blocco di A
                    &B[start_j],                 // intera matrice B
                    &C[start_i * ldC + start_j], // Offset del blocco di C
                    block_rows, N2, block_cols,
                    dbA, dbB, dbC);
    }
}

void matmatdist(MPI_Comm Gridcom, int ldA, int ldB, int ldC,
                double *A, double *B, double *C,
                int N1, int N2, int N3, int NPRow, int NPCol,
                int DB1, int DB2, int DB3, int NTrow, int NTcol)
{

    int rank, coords[2], size;
    MPI_Comm_size(Gridcom, &size);
    MPI_Comm_rank(Gridcom, &rank);
    MPI_Cart_coords(Gridcom, rank, 2, coords);

    int row = coords[0]; // Riga del processo
    int col = coords[1]; // Colonna del processo

    // Dimensioni dei blocchi locali
    int local_N1 = N1 / NPRow;
    int local_N2 = N2 / NPCol;
    int local_N3 = N3 / NPCol;

    // Allocazione dei blocchi locali
    double *local_A = (double *)malloc(local_N1 * local_N2 * sizeof(double));
    double *local_B = (double *)malloc(local_N2 * local_N3 * sizeof(double));
    double *local_C = (double *)calloc(local_N1 * local_N3, sizeof(double));

    // Scatter iniziale di A e B ai processi
    MPI_Scatter(A, local_N1 * local_N2, MPI_DOUBLE, local_A, local_N1 * local_N2, MPI_DOUBLE, 0, Gridcom);
    MPI_Scatter(B, local_N2 * local_N3, MPI_DOUBLE, local_B, local_N2 * local_N3, MPI_DOUBLE, 0, Gridcom);

    // Buffer per il broadcast
    double *Acol = (double *)malloc(local_N1 * local_N2 * sizeof(double));
    double *Brow = (double *)malloc(local_N2 * local_N3 * sizeof(double));

    for (int k = 0; k < NPCol; k++)
    {
        // Calcolo delle coordinate per il broadcast
        int source_col = (col - k + NPCol) % NPCol; // Broadcast lungo la riga
        int source_row = (row - k + NPRow) % NPRow; // Broadcast lungo la colonna

        // Broadcast del blocco di A lungo la riga
        if (col == k)
        {
            memcpy(Acol, local_A, local_N1 * local_N2 * sizeof(double));
        }
        MPI_Bcast(Acol, local_N1 * local_N2, MPI_DOUBLE, k, Gridcom);

        // Broadcast del blocco di B lungo la colonna
        if (row == k)
        {
            memcpy(Brow, local_B, local_N2 * local_N3 * sizeof(double));
        }
        MPI_Bcast(Brow, local_N2 * local_N3, MPI_DOUBLE, k, Gridcom);

        // Calcolo del prodotto locale
        matmatthread(ldA, ldB, ldC, Acol, Brow, local_C,
                     local_N1, local_N2, local_N3,
                     DB1, DB2, DB3, NTrow, NTcol);
    }

    // Raccolta dei risultati finali
    MPI_Gather(local_C, local_N1 * local_N3, MPI_DOUBLE, C, local_N1 * local_N3, MPI_DOUBLE, 0, Gridcom);

    // Pulizia
    free(local_A);
    free(local_B);
    free(local_C);
    free(Acol);
    free(Brow);
}

double get_cur_time()
{
    struct timeval tv;
    struct timezone tz;
    double cur_time;

    gettimeofday(&tv, &tz);
    cur_time = tv.tv_sec + tv.tv_usec / 1000000.0;

    return cur_time;
}