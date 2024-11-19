# Parallel and Distributed Computing Exercises

This repository contains a collection of exercises and examples for practicing **parallel** and **distributed computing** using **OpenMP** and **MPI** in **C**.  
The exercises are designed to demonstrate various concepts and techniques for efficient parallelization and inter-process communication.

---

## Features

- **OpenMP** examples for multi-threaded parallelism on shared-memory architectures.
- **MPI** examples for message-passing parallelism on distributed-memory systems.
- Hands-on exercises with increasing levels of complexity:
  - Simple parallel loops and reductions.
  - Load balancing and thread scheduling.
  - Collective communications, such as broadcast, scatter, and gather.
  - Point-to-point communications using MPI_Send and MPI_Recv.
  - Hybrid parallelism combining OpenMP and MPI.
- Detailed comments and explanations in the source code for learning purposes.

---

## Prerequisites

To compile and run the exercises, ensure the following tools and libraries are installed on your system:

- **C Compiler**: GCC or equivalent with support for OpenMP.
- **MPI Implementation**: [OpenMPI](https://www.open-mpi.org/) or [MPICH](https://www.mpich.org/).

### Check installation:

1. Verify OpenMP support:
   ```bash
   gcc -fopenmp -o test test.c
