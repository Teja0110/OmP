#include "mpi.h"
#include <omp.h>
#include <cstdlib>
#include <iostream>

#ifndef DEBUG
#define DEBUG false
#endif

using namespace std;

const int CELL_MAX = 4 + 1;

void populateMatrix(int * matrix, int n);
int* getTranspose(int* A, int n);
void printResult(int *matrix, int numprocs, int n);
void printMatrix(int *matrix, int n);

int main (int argc, char ** argv) {
  const int n = atoi(argv[1]);  //getting the input for the matrix size
  int numprocs, rank, numthreads;
  int *matrixA;
  int *matrixB;
  int *matrixC;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  matrixA = new int[n*n];
  matrixB = new int[n*n];
  matrixC = new int[n*n];

  populateMatrix(matrixA, n);
  if (DEBUG && rank == 0) {
    printMatrix(matrixA, n);
    cout << endl;
  }
  populateMatrix(matrixB, n);
  if (DEBUG && rank == 0) {
    printMatrix(matrixB, n);
    cout << endl;
  }


  double start, end;
  int tile = n/numprocs;
  int tilesize = n*n;
  int buffersize = tile * n;
  int sum;

  int *a = new int[buffersize];
  int *b = new int[buffersize];
  int *c = new int[tilesize];

  start = MPI_Wtime();
  matrixB = getTranspose(matrixB, n);

  MPI_Scatter(matrixA, buffersize, MPI_INT, a, buffersize, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Scatter(matrixB, buffersize, MPI_INT, b, buffersize, MPI_INT, 0, MPI_COMM_WORLD);

  # pragma omp parallel for            //OpenMP parallel for
  for (int i = 0; i < tile; i++) {
    numthreads = omp_get_num_threads();
    for (int j = 0; j < tile; j++) {
      sum = 0.0;
      for (int k = 0; k < n; k++) {
        sum += a[i*tile+k] * b[j*tile+k];
      }
      c[i*n+j] = sum;
    }
  }

  MPI_Gather(c, tilesize, MPI_INT, matrixC, tilesize, MPI_INT, 0, MPI_COMM_WORLD);

  end = MPI_Wtime(); // take the time after the multiplication
  double elapsed_time = end - start;

  if (DEBUG && rank == 0) {
    printResult(matrixC, numprocs, n);
  }

  if (rank == 0) {
    cout << "Num processes: " << numprocs << endl;
    cout << "Num threads: " << numthreads << endl;
    cout  << "Time for Multiplication: " << elapsed_time << endl;
  }

  //delete the matrices
    if (rank == 0) {
      delete[] matrixA;
      delete[] matrixB;
      delete[] matrixC;
    }

  MPI_Finalize();

  return 0;
}

void printResult(int *matrix, int numprocs, int n) {
  int tile = n / numprocs;
  int tilesize = tile * tile;

  int col = 0;
  for (int i = 0; i < numprocs; i++) {
    for (int j = 0; j < tilesize; j++) {
      cout << " " << matrix[i * n + j] << " \t";
      col++;
      if (col == n) {
        cout << endl;
        col = 0;
      }
    }
  }
}

void printMatrix(int *matrix, int n) {
  int col = 0;
  for (int i = 0; i < n * n; i++) {
    cout << " " << matrix[i] << " \t";
    col++;
    if (col ==n) {
      cout << endl;
      col = 0;
    }
  }
}

void populateMatrix(int * matrix, int n){
  for(int i = 0; i < n; i++ ){
    for(int j = 0; j < n; j ++){

      matrix[i*n+j] = rand() % CELL_MAX;
    }
  }
}

int* getTranspose(int* A, int n){
  int* matrix = new int[n*n]; //declare a new matrix to store the transpose of matrix B
  for(int i=0; i<n; i++){
    for(int j=0; j<n; j++){
      matrix[i*n+j] = A[j*n+i];
    }
  }
  return matrix;
}
