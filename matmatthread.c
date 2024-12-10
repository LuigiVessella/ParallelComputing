#include <omp.h>

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

        // range del blocco
        start_i = IDi * (N1 / NTrow);
        end_i = (IDi + 1) * (N1 / NTrow);

        start_j = IDj * (N3 / NTcol);
        end_j = (IDj + 1) * (N3 / NTcol);

        // blocchi da calcolare
        block_rows = end_i - start_i;
        block_cols = end_j - start_j;

        //printf("processo %d%d calcola %d %d della mat\n", IDi, IDj, block_rows, block_cols);

        matmatblock(ldA, ldB, ldC,
                    &A[start_i * ldA],           // Offset riga del blocco di A
                    B,                           // intera matrice B
                    &C[start_i * ldC + start_j], // Offset del blocco di C
                    block_rows, N2, block_cols,
                    dbA, dbB, dbC);
    }
}