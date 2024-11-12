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
            MPI_Irecv(daprev, LD, MPI_FLOAT, myid - 1, 0, MPI_COMM_WORLD, &request2); //request2 per daprev
        }

        

        // Se id != NP-1, invia l'ultima riga al processo successivo e ricevi la prima riga da questo
        if (myid != nproc - 1)
        {
            MPI_Isend(A + (N / nproc - 1) * LD, LD, MPI_FLOAT, myid + 1, 0, MPI_COMM_WORLD, &request3);
            MPI_Irecv(danext, LD, MPI_FLOAT, myid + 1, 0, MPI_COMM_WORLD, &request4); //request4 per danext
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
	    
            MPI_Wait(&request2, &status); //attendo daprev
            for (j = 1; j < N - 1; j++)
            {
                Anew[0 * LD + j] = 0.25 * (daprev[j] + A[1 * LD + j] +
                                           A[0 * LD + (j - 1)] + A[0 * LD + (j + 1)]);
            }
	
            MPI_Wait(&request1, &status);
       }




        // Calcola l'ultima riga di B se id != NP-1
        if (myid != nproc - 1)
        {

            MPI_Wait(&request4, &status); //attendo danext
            for (j = 1; j < N - 1; j++)
            {
                Anew[(N / nproc - 1) * LD + j] = 0.25 * (A[(N / nproc - 2) * LD + j] + danext[j] +
                                                         A[(N / nproc - 1) * LD + (j - 1)] + A[(N / nproc - 1) * LD + (j + 1)]);
            }

            MPI_Wait(&request3, &status);
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
