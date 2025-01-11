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

---

## Exercises

### Exercise 1: Parallel Sum with OpenMP
- **Description**: Implement a parallel sum of an array using OpenMP.
- **Concepts**: Parallel loops, reductions.
- **File**: `parallel_sum.c`

### Exercise 2: Matrix Multiplication with OpenMP
- **Description**: Implement matrix multiplication using OpenMP.
- **Concepts**: Parallel loops, nested parallelism.
- **File**: `matrix_multiplication.c`

### Exercise 3: Load Balancing with OpenMP
- **Description**: Explore different scheduling strategies in OpenMP.
- **Concepts**: Static, dynamic, and guided scheduling.
- **File**: `load_balancing.c`

### Exercise 4: MPI Hello World
- **Description**: Basic MPI program to print "Hello World" from multiple processes.
- **Concepts**: MPI initialization, rank, and size.
- **File**: `mpi_hello_world.c`

### Exercise 5: MPI Send and Receive
- **Description**: Implement point-to-point communication using MPI_Send and MPI_Recv.
- **Concepts**: Blocking communication, message passing.
- **File**: `mpi_send_receive.c`

### Exercise 6: MPI Collective Communication
- **Description**: Implement collective communication operations such as broadcast, scatter, and gather.
- **Concepts**: MPI_Bcast, MPI_Scatter, MPI_Gather.
- **File**: `mpi_collective.c`

### Exercise 7: Hybrid Parallelism with OpenMP and MPI
- **Description**: Combine OpenMP and MPI for hybrid parallelism.
- **Concepts**: Shared-memory and distributed-memory parallelism.
- **File**: `hybrid_parallelism.c`

### Exercise 8: Parallel Max Sum with OpenMP
- **Description**: Implement a parallel algorithm to find the maximum sum of elements in an array using OpenMP.
- **Concepts**: Parallel loops, critical sections.
- **File**: `maxsum.c`

---

## Compilation and Execution

### OpenMP Exercises
To compile an OpenMP exercise, use the following command:
```sh
gcc -fopenmp -o output_file source_file.c