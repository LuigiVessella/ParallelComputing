/*
con homebrew installare gcc (non necessariamente) e libomp e poi compilare con
clang -Xclang -fopenmp -L/opt/homebrew/opt/libomp/lib -I/opt/homebrew/opt/libomp/include -lomp maxsum.c -o maxsum
Mac Silicon M1
*/


// VESSELLA LUIGI N97000503
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <sys/time.h>

double maxsum(int, int, double *, int);

double get_cur_time();

int main()
{

    int NT, N, i, j, LD;
    double MAX = 0, *A;
    double t1, t2, save;

    LD = 800;
    A = (double *)malloc(sizeof(double) * LD * LD);
    N = 800;
    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            A[i * LD + j] = (rand() % 100);
        }
    }

    for (NT = 1; NT <= 8; NT = NT * 2)
    {

        printf("finalmente posso compilare OMP\n");
        printf("===============\n");
        printf(" NUMERO THREAD = %d \n", NT);

        t1 = get_cur_time();
        MAX = maxsum(N, LD, A, NT);
        t2 = get_cur_time();

        if (NT == 1)
            save = t2 - t1;
        printf("il massimo della somma dei moduli con N = %d e' %f \n", N, MAX);
        printf("il tempo totale e' %e , lo speedup = %f , l'efficienza = %f \n", t2 - t1, save / (t2 - t1), save / (t2 - t1) / NT);
    }
}

double maxsum(int N, int LD, double *A, int NT)
{
    double MAX = 0.0, temp_sum = 0.0;
    int i, j, id, start = 0, end = 0;

    omp_set_num_threads(NT);

    #pragma omp parallel private(i, j, id, temp_sum, start, end)
    {

        id = omp_get_thread_num();
        start = id * N / NT;
        end = (id + 1) * N / NT;

        for (i = start; i < end; i++)
        {
            temp_sum = 0;
            for (j = 0; j < N; j++)
            {
                temp_sum += sqrt(A[i * LD + j]);
            }

            #pragma omp critical
            if (temp_sum > MAX)
            {
                MAX = temp_sum;
            }
        }
    }
    return MAX;
}
double get_cur_time()
{
    struct timeval tv;
    struct timezone  tz;
    double cur_time;

    gettimeofday(&tv, &tz);
    cur_time = tv.tv_sec + tv.tv_usec / 1000000.0;

    return cur_time;
}