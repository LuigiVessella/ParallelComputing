#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

void matmatijk(int ldA, int ldB, int ldC, double *A, double *B, double *C, int N1, int N2, int N3);
void matmatjik(int ldA, int ldB, int ldC, double *A, double *B, double *C, int N1, int N2, int N3);
void matmatikj(int ldA, int ldB, int ldC, double *A, double *B, double *C, int N1, int N2, int N3);
void matmatjki(int ldA, int ldB, int ldC, double *A, double *B, double *C, int N1, int N2, int N3);
void matmatkij(int ldA, int ldB, int ldC, double *A, double *B, double *C, int N1, int N2, int N3);
void matmatkji(int ldA, int ldB, int ldC, double *A, double *B, double *C, int N1, int N2, int N3);
void matmatblock(int ldA, int ldB, int ldC, double *A, double *B, double *C, int N1, int N2, int N3, int dbA, int dbB, int dbC);

double get_cur_time();

int main()
{
    int N1 = 512, N2 = 512, N3 = 512;
    int ldA = 1000, ldB = 1000, ldC = 1000;
    double *A, *B, *C;
    int i, j, k;
    int best = 0;
    double gflops = 0, gflops2 = 0;

    // A matrice di N1 righe e N2 colonne
    // B matrice di N2 righe e N3 colonne
    // C matrice di N1 righe e N3 colonne

    A = (double *)malloc(1000 * 1000 * sizeof(double));
    B = (double *)malloc(1000 * 1000 * sizeof(double));
    C = (double *)malloc(1000 * 1000 * sizeof(double));

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

    double t1 = get_cur_time();
    matmatijk(ldA, ldB, ldC, A, B, C, N1, N2, N3);
    double t2 = get_cur_time();
    gflops = (2 * (N1 * N2 * N3)) / (t2 - t1) / 1e9;
    printf("time ijk %f gflops ijk %f\n", t2 - t1, gflops);

    best = gflops;

    t1 = get_cur_time();
    matmatjik(ldA, ldB, ldC, A, B, C, N1, N2, N3);
    t2 = get_cur_time();
    gflops2 = (2 * (N1 * N2 * N3)) / (t2 - t1) / 1e9;
    printf("time jik %f gflops jik %f\n", t2 - t1, gflops2);

    if (gflops2 < best)
        best = gflops2;

    t1 = get_cur_time();
    matmatikj(ldA, ldB, ldC, A, B, C, N1, N2, N3);
    t2 = get_cur_time();
    gflops2 = (2 * (N1 * N2 * N3)) / (t2 - t1) / 1e9;
    printf("time ikj %f gflops ikj %f\n", t2 - t1, gflops2);

    if (gflops2 < best)
        best = gflops2;

    t1 = get_cur_time();
    matmatjki(ldA, ldB, ldC, A, B, C, N1, N2, N3);
    t2 = get_cur_time();
    gflops2 = (2 * (N1 * N2 * N3)) / (t2 - t1) / 1e9;
    printf("time jki %f gflops jki %f\n", t2 - t1, gflops2);

    if (gflops2 < best)
        best = gflops2;

    t1 = get_cur_time();
    matmatkij(ldA, ldB, ldC, A, B, C, N1, N2, N3);
    t2 = get_cur_time();
    gflops2 = (2 * (N1 * N2 * N3)) / (t2 - t1) / 1e9;
    printf("time kij %f gflops kij %f\n", t2 - t1, gflops2);

    if (gflops2 < best)
        best = gflops2;

    t1 = get_cur_time();
    matmatkji(ldA, ldB, ldC, A, B, C, N1, N2, N3);
    t2 = get_cur_time();
    gflops2 = (2 * (N1 * N2 * N3)) / (t2 - t1) / 1e9;
    printf("time kji %f gflops kji %f\n", t2 - t1, gflops2);

    if (gflops2 < best)
        best = gflops2;

    t1 = get_cur_time();
    matmatblock(ldA, ldB, ldC, A, B, C, N1, N2, N3, 512, 512, 512);
    t2 = get_cur_time();
    gflops2 = (2 * (N1 * N2 * N3)) / (t2 - t1) / 1e9;
    printf("time matmatblock kji %f gflops kji %f\n", t2 - t1, gflops2);

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

void matmatblock(int ldA, int ldB, int ldC, double *A, double *B, double *C, int N1, int N2, int N3, int dbA, int dbB, int dbC)
{
    int i, j, k;
    int ii, jj, kk;

    for (ii = 0; ii < N1; ii += dbA)
    {
        for (kk = 0; kk < N3; kk += dbC)
        {
            for (jj = 0; jj < N2; jj += dbB)
            {
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

double get_cur_time()
{
    struct timeval tv;
    struct timezone tz;
    double cur_time;

    gettimeofday(&tv, &tz);
    cur_time = tv.tv_sec + tv.tv_usec / 1000000.0;

    return cur_time;
}