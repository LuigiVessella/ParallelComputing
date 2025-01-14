/*
Per compilare, su Mac silicon M1, è necessario utilizzare il comando:
clang -O3 -Xclang -fopenmp -L/opt/homebrew/opt/libomp/lib -I/opt/homebrew/opt/libomp/include -lomp matmathread.c -o matmathread

Questo dopo aver installato libomp con il comando:
brew install libomp

Su altre piattaforme, la procedura è simile, ma potrebbero essere necessarie modifiche ai comandi.

Ricordare la flag -O3 per l'ottimizzazione del codice, altrimenti non otterrete i gflops attesi.
*/

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <string.h>
#include <omp.h>

void matmatijk(int ldA, int ldB, int ldC, double *A, double *B, double *C, int N1, int N2, int N3);
void matmatjik(int ldA, int ldB, int ldC, double *A, double *B, double *C, int N1, int N2, int N3);
void matmatikj(int ldA, int ldB, int ldC, double *A, double *B, double *C, int N1, int N2, int N3);
void matmatjki(int ldA, int ldB, int ldC, double *A, double *B, double *C, int N1, int N2, int N3);
void matmatkij(int ldA, int ldB, int ldC, double *A, double *B, double *C, int N1, int N2, int N3);
void matmatkji(int ldA, int ldB, int ldC, double *A, double *B, double *C, int N1, int N2, int N3);
void initialize_matrices(int ldA, int ldB, int ldC, double *A, double *B, double *C, int N1, int N2, int N3);
double calculate_gflops(int N1, int N2, int N3, double time_seconds);
int compare_matrices(double *C1, double *C2, int ldC, int N1, int N3);
void matmatblock(int ldA, int ldB, int ldC, double *A, double *B, double *C, int N1, int N2, int N3, int dbA, int dbB, int dbC);
void matmatthread(int ldA, int ldB, int ldC, double *A, double *B, double *C, int N1, int N2, int N3, int dbA, int dbB, int dbC, int NTROW, int NTCOL);
double get_cur_time();

int main()
{
    int dimensions[] = {1024, 2048};
    int num_dimensions = 2;
    int ldA = 3500, ldB = 3500, ldC = 3500;
    int dbA = 256, dbB = 256, dbC = 256;
    double t1, t2;
    double *A, *B, *C, *C_block;
    int i, j, k, dim_idx;
    double gflops = 0.0, gflops2 = 0.0;
    double cicler = 0.0;

    int NTROW_values[] = {1, 2, 2, 2};
    int NTCOL_values[] = {1, 1, 2, 4};
    int num_combinations = 4;
    int idx;

    for (dim_idx = 0; dim_idx < num_dimensions; dim_idx++)
    {
        int N1 = dimensions[dim_idx];
        int N2 = dimensions[dim_idx];
        int N3 = dimensions[dim_idx];

        printf("\n\n=== Test con dimensione %dx%dx%d ===\n", N1, N2, N3);

        C_block = (double *)malloc(ldC * ldC * sizeof(double));
        A = (double *)malloc(ldA * ldA * sizeof(double));
        B = (double *)malloc(ldB * ldB * sizeof(double));
        C = (double *)malloc(ldC * ldC * sizeof(double));

        for (idx = 0; idx < num_combinations; idx++)
        {
            int NTROW = NTROW_values[idx];
            int NTCOL = NTCOL_values[idx];

            printf("\n----------- Iterazione %d: NTROW = %d, NTCOL = %d -----------\n", idx + 1, NTROW, NTCOL);

            // inizializzo matrici
            initialize_matrices(ldA, ldB, ldC, A, B, C, N1, N2, N3);

            t1 = get_cur_time();
            matmatijk(ldA, ldB, ldC, A, B, C, N1, N2, N3);
            t2 = get_cur_time();
            gflops = calculate_gflops(N1, N2, N3, t2 - t1);
            printf("ijk: time = %f s, GFLOPS = %f\n", t2 - t1, gflops);

            // Copia il risultato di C in C_block per il confronto
            memcpy(C_block, C, ldC * ldC * sizeof(double));

            // Esegui matmatblock
            initialize_matrices(ldA, ldB, ldC, A, B, C, N1, N2, N3);
            t1 = get_cur_time();
            matmatblock(ldA, ldB, ldC, A, B, C, N1, N2, N3, dbA, dbB, dbC);
            t2 = get_cur_time();
            gflops2 = calculate_gflops(N1, N2, N3, t2 - t1);
            printf("block: time = %f s, GFLOPS = %f\n", t2 - t1, gflops2);

            int check = compare_matrices(C, C_block, ldC, N1, N3);

            // Confronta i risultati
            if (check)
            {
                printf("Le matrici sono uguali.\n");
            }
            else
            {
                printf("Le matrici non sono uguali.\n");
            }

            // Esegui matmatthread
            initialize_matrices(ldA, ldB, ldC, A, B, C, N1, N2, N3);
            t1 = get_cur_time();
            matmatthread(ldA, ldB, ldC, A, B, C, N1, N2, N3, dbA, dbB, dbC, NTROW, NTCOL);
            t2 = get_cur_time();
            gflops2 = calculate_gflops(N1, N2, N3, t2 - t1);
            printf("thread (NTROW=%d, NTCOL=%d): time = %f s, GFLOPS = %f\n", NTROW, NTCOL, t2 - t1, gflops2);
            check = compare_matrices(C, C_block, ldC, N1, N3);
            // Confronta i risultati
            if (check)
            {
                printf("Le matrici sono uguali.\n");
            }
            else
            {
                printf("Le matrici non sono uguali.\n");
            }
        }
        free(C_block);
        free(A);
        free(B);
        free(C);
    }

    return 0;
    /*
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
        */

    /*
        // versione parallela
        printf("versione parallela 1-1\n");
        t1 = get_cur_time();
        matmatthread(ldA, ldB, ldC, A, B, C, N1, N2, N3, 2, 2, 2, 1, 1);
        t2 = get_cur_time();
        gflops2 = ((2.0 * ((double)N1 * (double)N2 * (double)N3)) / (t2 - t1)) / 1e9;
        printf("time versione parallela 1-1 %f gflops ikj %f\n", t2 - t1, gflops2);

        // versione parallela
        printf("versione parallela 1-2\n");
        t1 = get_cur_time();
        matmatthread(ldA, ldB, ldC, A, B, C, N1, N2, N3, 2, 2, 2, 1, 2);
        t2 = get_cur_time();
        gflops2 = ((2.0 * ((double)N1 * (double)N2 * (double)N3)) / (t2 - t1)) / 1e9;
        printf("time versione parallela 1-2 %f gflops ikj %f\n", t2 - t1, gflops2);
    */

    /*
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
    */
    // versione parallela
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

        // printf("processo %d%d calcola %d %d della mat\n", IDi, IDj, block_rows, block_cols);
#pragma omp critical
        {
            printf("\n=== Thread %d (IDi=%d, IDj=%d) ===\n", thread_id, IDi, IDj);
            printf("Divisione teorica: N1/NTrow=%d, N3/NTcol=%d\n", N1 / NTrow, N3 / NTcol);
            printf("Range righe: start_i=%d, end_i=%d (block_rows=%d)\n", start_i, end_i, block_rows);
            printf("Range colonne: start_j=%d, end_j=%d (block_cols=%d)\n", start_j, end_j, block_cols);
        }

        // Chiamata alla funzione `matmatblock` per calcolare il blocco
        matmatblock(ldA, ldB, ldC,
                    &A[start_i * ldA],           // Offset riga del blocco di A
                    &B[start_j],                 // Mi serve l'intera colonna di B
                    &C[start_i * ldC + start_j], // Offset del blocco di C
                    block_rows, N2, block_cols,
                    dbA, dbB, dbC);
    }
}

void initialize_matrices(int ldA, int ldB, int ldC, double *A, double *B, double *C, int N1, int N2, int N3)
{
    int i, j;
    for (i = 0; i < N1; i++)
    {
        for (j = 0; j < N2; j++)
        {
            A[i * ldA + j] = 2;
        }
    }

    for (i = 0; i < N2; i++)
    {
        for (j = 0; j < N3; j++)
        {
            B[i * ldB + j] = 3;
        }
    }

    for (i = 0; i < N1; i++)
    {
        for (j = 0; j < N3; j++)
        {
            C[i * ldC + j] = 0;
        }
    }
}

int compare_matrices(double *C1, double *C2, int ldC, int N1, int N3)
{
    int i, j;
    for (i = 0; i < N1; i++)
    {
        for (j = 0; j < N3; j++)
        {
            if (C1[i * ldC + j] != C2[i * ldC + j])
            {
                return 0; // le matrici non sono uguali
            }
        }
    }
    return 1; // le matrici sono uguali
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
double calculate_gflops(int N1, int N2, int N3, double time_seconds)
{
    double operations = (double)(2.0 * ((double)N1 * (double)N2 * (double)N3));
    double gflops = (double)((operations / time_seconds) / 1e9);
    return gflops;
}