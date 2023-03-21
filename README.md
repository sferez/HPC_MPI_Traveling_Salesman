
# High Performance Technical Computing Project

This project focuses on the implementation of a High-Performance Technical Computing solution on an HPC cluster, utilizing MPI to parallelize the wandering salesman problem.

## Report

To be added.

## Structure

The following code folders are included in this project:

- 'BruteForce' contains the brute force implementation of the Wandering Salesman Problem
- 'Branch_and_Bound' contains the branch and bound implementation of the Wandering Salesman Problem
- 'input' contains the input files for the algorithms from 4 to 18 cities

Both folders contain a serial and a parallel implementation of the algorithm.

# Usage

- Have MPI installed on your system.
- Compile the code using mpic++.
- Run the code using mpirun.

## Example

### Installing MPI
```bash
brew install open-mpi # For MacOS
```

### Compiling
```bash
mpic++ Branch_and_Bound/mainMPI.cpp # Compiles the code
# or
mpic++ BruteForce/mainMPI.cpp # Compiles the code
```
Get the executable file 'a.out'.

### Running
```bash
# Run the code on 4 processes
# using the input file 'input/dist16'
mpirun -np 4 a.out input/dist16 
# or
# Run the code on 8 processes
# using the input file 'input/dist17'
mpirun -np 8 a.out input/dist17
```

## Authors

- [@sferez](https://github.com/sferez)
