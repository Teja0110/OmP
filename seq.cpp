#include <iostream>
#include <cstdlib> 
#include <ctime>
#include <cmath>
#include <omp.h>

using namespace std;

void populateMatrix(double ** matrix, int n);
double ** initializeMatrix(int n);
double multiplyMatrix(int n);
double** getTranspose(double** A, int n);
 
void populateMatrix(double ** matrix, int n){

	for(int i = 0; i < n; i++ ){
		for(int j = 0; j < n; j ++){
			matrix[i][j] = rand();
		}
	}
}

 
double ** initializeMatrix(int n){
	double ** matrix;
	matrix = new double*[n];;
	for (int row = 0; row<n; row++) {
	  matrix[row] = new double[n];
	}
	return matrix;
}
 
double** getTranspose(double** A, int n){
	double** matrix = new double*[n]; //declare a new matrix to store the transpose of matrix B
	for(int i=0; i<n; i++){  
		matrix[i]=new double[n];
		for(int j=0; j<n; j++){
			matrix[i][j] = A[j][i];
		}
	}
	return matrix;
}

double multiplyMatrix(int n) {

 
	double **matrixA;
	double **matrixB;
	double **matrixC;
	double **matrixBT;

 
	matrixA = initializeMatrix(n);
	matrixB = initializeMatrix(n);
	matrixC = initializeMatrix(n);
	matrixBT = initializeMatrix(n);
 
	populateMatrix(matrixA, n);
	populateMatrix(matrixB, n);
        
    double start, end;
	start = omp_get_wtime();   
	/*
	*matrix multiplication (parallel_optimized)
	*/
	matrixBT = getTranspose(matrixB, n);  //get the transpose of matrix B
	double sum;
	cout << "Threads 1: " << omp_get_num_threads() << endl;
	
 	//OpenMP parallel for
	for(int i = 0; i < n; i++){
	        for(int j = 0; j < n; j++){
	            sum = 0.0;
	            for (int k = 0; k < n; k++) {
	            	sum += matrixA[i][k] *  matrixBT[j][k];;
	            }
	            matrixC[i][j] =sum;
	        }
	    }
	/*
	*end of matrix multiplication (parallel_optimized)
	*/

	end = omp_get_wtime(); //take the time after the multiplication
	double elapsed_time = end -start;

	//delete the matrices
	for (size_t i = 0; i < n; i++) {
        delete [] matrixA[i];
        delete [] matrixB[i];
        delete [] matrixC[i];
    }
	delete[] matrixA;
	delete[] matrixB;
	delete[] matrixC;

	return elapsed_time;
}
 
int main(int argc, const char * argv[]) {
	int n = atoi(argv[1]);  //getting the input for the matrix size
	double elapsed_time = multiplyMatrix(n);
	cout  << "Time for Multiplication with matrix transpose: " << elapsed_time << endl;
 
	return 0;
}