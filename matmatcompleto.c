#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <omp.h>

void matmatijk(int ldA, int ldB, int ldC, double *A, double *B, double *C, int N1, int N2, int N3);
void matmatjik(int ldA, int ldB, int ldC, double *A, double *B, double *C, int N1, int N2, int N3);
void matmatikj(int ldA, int ldB, int ldC, double *A, double *B, double *C, int N1, int N2, int N3);
void matmatjki(int ldA, int ldB, int ldC, double *A, double *B, double *C, int N1, int N2, int N3);
void matmatkij(int ldA, int ldB, int ldC, double *A, double *B, double *C, int N1, int N2, int N3);
void matmatkji(int ldA, int ldB, int ldC, double *A, double *B, double *C, int N1, int N2, int N3);

void matmatblock(int ldA, int ldB, int ldC, double *A, double *B, double *C, int N1, int N2, int N3, int dbA, int dbB, int dbC);

void matmatthread(int ldA, int ldB, int ldC, double *A, double *B, double *C, int N1, int N2, int N3, int dbA, int dbB, int dbC, int NTROW, int NTCOL);
double get_cur_time();

int main()
{

    int N1 = 8, N2 = 8, N3 = 8;
    int ldA = 50, ldB = 50, ldC = 50;
    double t1, t2;
    double *A, *B, *C;
    int i, j, k;
    int best = 0;
    double gflops = 0.0, gflops2 = 0.0;
    double cicler = 0.0;

    // A matrice di N1 righe e N2 colonne
    // B matrice di N2 righe e N3 colonne
    // C matrice di N1 righe e N3 colonne

    A = (double *)malloc(ldA * ldA * sizeof(double));
    B = (double *)malloc(ldB * ldB * sizeof(double));
    C = (double *)malloc(ldC * ldC * sizeof(double));

    for (i = 0; i < N1; i++)
    {
        for (j = 0; j < N2; j++)
        {
            // usiamo lo stesso ciclo in quanto le dimensioni sono assolutamente uguali
            A[i * ldA + j] = 2;
            B[i * ldB + j] = 3;
            C[i * ldC + j] = 0;
        }
    }

    printf("B:\n\n");
    for (i = 0; i < N1; i++)
    {
        for (j = 0; j < N2; j++)
        {
            // usiamo lo stesso ciclo in quanto le dimensioni sono assolutamente uguali
            printf("%f ", B[i * ldC + j]);
        }
        printf("\n");
    }

    matmatijk(ldA, ldB, ldC, A, B, C, N1, N2, N3);

    for (i = 0; i < N1; i++)
    {
        for (j = 0; j < N2; j++)
        {
            // usiamo lo stesso ciclo in quanto le dimensioni sono assolutamente uguali
            printf("%f ", C[i * ldC + j]);
        }
        printf("\n");
    }

    for (i = 0; i < N1; i++)
    {
        for (j = 0; j < N2; j++)
        {
            // usiamo lo stesso ciclo in quanto le dimensioni sono assolutamente uguali
            A[i * ldA + j] = 2;
            B[i * ldB + j] = 3;
            C[i * ldC + j] = 0;
        }
    }

    printf("\ns-------------------------\n\n\n");
    matmatblock(ldA, ldB, ldC, A, B, C, N1, N2, N3, 2, 2, 2);

    for (i = 0; i < N1; i++)
    {
        for (j = 0; j < N2; j++)
        {
            // usiamo lo stesso ciclo in quanto le dimensioni sono assolutamente uguali
            printf("%f ", C[i * ldC + j]);
        }
        printf("\n");
    }

    for (i = 0; i < N1; i++)
    {
        for (j = 0; j < N2; j++)
        {
            // usiamo lo stesso ciclo in quanto le dimensioni sono assolutamente uguali
            A[i * ldA + j] = 2;
            B[i * ldB + j] = 3;
            C[i * ldC + j] = 0;
        }
    }

    printf("\ns------------MATMATTHREAD-------------\n\n\n");
    matmatthread(ldA, ldB, ldC, A, B, C, N1, N2, N3, 2, 2, 4, 2, 2);

    for (i = 0; i < N1; i++)
    {
        for (j = 0; j < N2; j++)
        {
            // usiamo lo stesso ciclo in quanto le dimensioni sono assolutamente uguali
            printf("%f ", C[i * ldC + j]);
        }
        printf("\n");
    }



        // inizializzo matrici
        for (i = 0; i < N1; i++)
        {
            for (j = 0; j < N2; j++)
            {
                // usiamo lo stesso ciclo in quanto le dimensioni sono assolutamente uguali
                A[i * ldA + j] = j;
                B[i * ldB + j] = j;
                C[i * ldC + j] = 0;
            }
        }

        t1 = get_cur_time();
        matmatijk(ldA, ldB, ldC, A, B, C, N1, N2, N3);
        t2 = get_cur_time();
        gflops = ((2.0 * ((double)N1 * (double)N2 * (double)N3)) / (t2 - t1)) / 1e9;
        printf("time ijk %f gflops ijk %f\n", t2 - t1, gflops);

        best = gflops;

        t1 = get_cur_time();
        matmatjik(ldA, ldB, ldC, A, B, C, N1, N2, N3);
        t2 = get_cur_time();
        gflops2 = ((2.0 * ((double)N1 * (double)N2 * (double)N3)) / (t2 - t1)) / 1e9;
        printf("time jik %f gflops jik %f\n", t2 - t1, gflops2);

        if (gflops2 < best)
            best = gflops2;

        t1 = get_cur_time();
        matmatikj(ldA, ldB, ldC, A, B, C, N1, N2, N3);
        t2 = get_cur_time();
        gflops2 = ((2.0 * ((double)N1 * (double)N2 * (double)N3)) / (t2 - t1)) / 1e9;
        printf("time ikj %f gflops ikj %f\n", t2 - t1, gflops2);

        if (gflops2 < best)
            best = gflops2;

        t1 = get_cur_time();
        matmatjki(ldA, ldB, ldC, A, B, C, N1, N2, N3);
        t2 = get_cur_time();
        gflops2 = ((2.0 * ((double)N1 * (double)N2 * (double)N3)) / (t2 - t1)) / 1e9;
        printf("time jki %f gflops jki %f\n", t2 - t1, gflops2);

        if (gflops2 < best)
            best = gflops2;

        t1 = get_cur_time();
        matmatkij(ldA, ldB, ldC, A, B, C, N1, N2, N3);
        t2 = get_cur_time();
        gflops2 = ((2.0 * ((double)N1 * (double)N2 * (double)N3)) / (t2 - t1)) / 1e9;
        printf("time kij %f gflops kij %f\n", t2 - t1, gflops2);

        if (gflops2 < best)
            best = gflops2;

        t1 = get_cur_time();
        matmatkji(ldA, ldB, ldC, A, B, C, N1, N2, N3);
        t2 = get_cur_time();
        gflops2 = ((2.0 * ((double)N1 * (double)N2 * (double)N3)) / (t2 - t1)) / 1e9;
        printf("time kji %f gflops kji %f\n", t2 - t1, gflops2);

        if (gflops2 < best)
            best = gflops2;

        t1 = get_cur_time();
        matmatblock(ldA, ldB, ldC, A, B, C, N1, N2, N3, 64, 64, 64);
        t2 = get_cur_time();
        gflops2 = ((2.0 * ((double)N1 * (double)N2 * (double)N3)) / (t2 - t1)) / 1e9;
        printf("time matmatblock ikj %f gflops ikj %f\n", t2 - t1, gflops2);



        // versione parallela
        printf("versione parallela\n\n");
        t1 = get_cur_time();
        matmatthread(ldA, ldB, ldC, A, B, C, N1, N2, N3, 64, 64, 64, 2, 2);
        t2 = get_cur_time();
        gflops2 = ((2.0 * ((double)N1 * (double)N2 * (double)N3)) / (t2 - t1)) / 1e9;
        printf("time matmatblockv2 ikj %f gflops ikj %f\n", t2 - t1, gflops2);
    }

    return 0;
}

void matmatijk(int ldA, int ldB, int ldC, double *A, double *B, double *C, int N1, int N2, int N3)
{

    int i, j, k;
    for (i = 0; i < N1; i++)
    {
        for (j = 0; j < N2; j++)
        {
            for (k = 0; k < N3; k++)
            {
                C[i * ldC + j] = C[i * ldC + j] + (A[i * ldA + k] * B[k * ldB + j]);
            }
        }
    }
}

void matmatjik(int ldA, int ldB, int ldC, double *A, double *B, double *C, int N1, int N2, int N3)
{

    int i, j, k;
    for (j = 0; j < N1; j++)
    {
        for (i = 0; i < N2; i++)
        {
            for (k = 0; k < N3; k++)
            {
                C[i * ldC + j] = C[i * ldC + j] + (A[i * ldA + k] * B[k * ldB + j]);
            }
        }
    }
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

void matmatjki(int ldA, int ldB, int ldC, double *A, double *B, double *C, int N1, int N2, int N3)
{

    int i, j, k;
    for (j = 0; j < N1; j++)
    {
        for (k = 0; k < N2; k++)
        {
            for (i = 0; i < N3; i++)
            {
                C[i * ldC + j] = C[i * ldC + j] + (A[i * ldA + k] * B[k * ldB + j]);
            }
        }
    }
}

void matmatkij(int ldA, int ldB, int ldC, double *A, double *B, double *C, int N1, int N2, int N3)
{

    int i, j, k;
    for (k = 0; k < N1; k++)
    {
        for (i = 0; i < N2; i++)
        {
            for (j = 0; j < N3; j++)
            {
                C[i * ldC + j] = C[i * ldC + j] + (A[i * ldA + k] * B[k * ldB + j]);
            }
        }
    }
}

void matmatkji(int ldA, int ldB, int ldC, double *A, double *B, double *C, int N1, int N2, int N3)
{

    int i, j, k;
    for (k = 0; k < N1; k++)
    {
        for (j = 0; j < N2; j++)
        {
            for (i = 0; i < N3; i++)
            {
                C[i * ldC + j] = C[i * ldC + j] + (A[i * ldA + k] * B[k * ldB + j]);
            }
        }
    }
}
/*
void matmatblockv2(int ldA, int ldB, int ldC, double *A, double *B, double *C, int N1, int N2, int N3, int dbA, int dbB, int dbC)
{
    int i, j, k;
    int ii, jj, kk;

    for (ii = 0; ii < N1; ii += dbA)
    {
        for (kk = 0; kk < N3; kk += dbC)
        {
            for (jj = 0; jj < N2; jj += dbB)
            {
                // matmatikj(ldA, ldB, ldC, &A[ii * ldA + kk], &B[kk * ldB + jj], &C[ii * ldC + jj], ii + dbA, jj + dbB, kk + dbC );
                // Cicli interni per calcolare i blocchi
                for (i = ii; i < ii + dbA && i < N1; i++)
                {
                    for (k = kk; k < kk + dbC && k < N3; k++)
                    {
                        for (j = jj; j < jj + dbB && j < N2; j++)
                        {
                            C[i * ldC + j] += A[i * ldA + k] * B[k * ldB + j];
                        }
                    }
                }
            }
        }
    }
}
*/

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
                    B,                           // intera matrice B
                    &C[start_i * ldC + start_j], // Offset del blocco di C
                    block_rows, N2, block_cols,
                    dbA, dbB, dbC);
    }
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