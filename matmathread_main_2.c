/*
Come il file matmathread.c ma con un main differente, i risultati prodotti sono gli stessi.
*/

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <math.h>
#include <omp.h>


void initialize_matrix_to_zero(double *A, int LD, int N);
void initialize_matrix(double *A, int LD, int N);
int compare_matrix(double *A, double *B, int N, int ld);
void matmatblock(int ldA, int ldB, int ldC, double *A, double *B, double *C, int N1, int N2, int N3, int dbA, int dbB, int dbC);
void matmatthread(int ldA, int ldB, int ldC, double *A, double *B, double *C, int N1, int N2, int N3, int dbA, int dbB, int dbC, int NTROW, int NTCOL);
void matmatikj(int ldA, int ldB, int ldC, double *A, double *B, double *C, int N1, int N2, int N3);

double get_cur_time();

int main(int argc, char *argv[]) {

  double *A, *B, *C, *C_ref; 
  int N1, N2, N3, ldA, ldB, ldC, dbA, dbB, dbC, NTrow, NTcol;
  int i, j;
  double t1, t2, t_seq, t_par;
  double operations_number;
  double gflops_matmatblock, gflops_matmatthread, speedup, efficiency;
  int combinations[][2] = {{1,1}, {1,2}, {2,2}, {2,4}};

  ldA = 3500;
  ldB = 3500;
  ldC = 3500;

  A = (double*)malloc(ldA * ldA * sizeof(double));
  B = (double*)malloc(ldB * ldB * sizeof(double));
  C = (double*)malloc(ldC * ldC * sizeof(double));
  C_ref = (double*)malloc(ldC * ldC * sizeof(double));

  dbA = 256;
  dbB = 256;
  dbC = 256;

  for (i = 0; i < 2; i++) {
    N1 = pow(2, 10 + i);
    N2 = pow(2, 10 + i);
    N3 = pow(2, 10 + i);

    initialize_matrix(A, ldA, N1);
    initialize_matrix(B, ldB, N2);

    operations_number = 2.00 * (double)N1 * (double)N2 * (double)N3;

    initialize_matrix_to_zero(C_ref, ldC, N3); 
    t1 = get_cur_time();
    matmatblock(ldA, ldB, ldC, A, B, C_ref, N1, N2, N3, dbA, dbB, dbC); 
    t2 = get_cur_time();
    t_seq = t2 - t1;
    gflops_matmatblock = (operations_number / t_seq) / 1e9;

    printf("-------------------------------\n");
    printf("Matrix size: %d x %d\n", N1, N1);
    printf("matmatblock: %f Gflops\n\n", gflops_matmatblock); 

    for (j = 0; j < 4; j++) {
      NTrow = combinations[j][0]; 
      NTcol = combinations[j][1];

      initialize_matrix_to_zero(C, ldC, N3);
      t1 = get_cur_time();
      matmatthread(ldA, ldB, ldC, A, B, C, N1, N2, N3, dbA, dbB, dbC, NTrow, NTcol);
      t2 = get_cur_time();
      t_par = t2 - t1;
      gflops_matmatthread = (operations_number / t_par) / 1e9;

      speedup = t_seq / t_par;
      efficiency = speedup / (NTrow * NTcol);

      if (!compare_matrix(C, C_ref, N3, ldC)){ 
        printf("NT = %d\n", NTrow * NTcol);
        printf("The matrices are different, exiting.\n\n"); 
        return 0;
      }

      printf("NT = %d\n", NTrow * NTcol); 
      printf("matmatthread: %f Gflops\n", gflops_matmatthread); 
      printf("Speedup: %f\n", speedup); 
      printf("Efficiency: %f\n\n", efficiency);
    }
  }

  free(A);
  free(B);
  free(C);
  free(C_ref);

  return 0;
}

void initialize_matrix(double *A, int LD, int N){
  int i, j;

  for(i = 0; i < N; i++){
    for(j = 0; j < N; j++){
      A[i * LD + j] = rand() % 100;
    }
  }
}

void initialize_matrix_to_zero(double *A, int LD, int N){
  int i, j;

  for(i = 0; i < N; i++){
    for(j = 0; j < N; j++){
      A[i * LD + j] = 0;
    }
  }
}

int compare_matrix(double *A, double *B, int N, int ld){
  int i, j;

  for(i = 0; i < N; i++){
    for(j = 0; j < N; j++){
      if(A[i * ld + j] != B[i * ld + j]){
        return 0;
      }
    }
  }

  return 1;
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

        // Chiamata alla funzione `matmatblock` per calcolare il blocco

        matmatblock(ldA, ldB, ldC,
                    &A[start_i * ldA],           // Offset riga del blocco di A
                    &B[start_j],                 // intera matrice B
                    &C[start_i * ldC + start_j], // Offset del blocco di C
                    block_rows, N2, block_cols,
                    dbA, dbB, dbC);
    }
}

double get_cur_time() {

  struct timeval   tv;
  struct timezone  tz;
  double cur_time;

  gettimeofday(&tv, &tz);
  cur_time = tv.tv_sec + tv.tv_usec / 1000000.0;

  return cur_time;
}
