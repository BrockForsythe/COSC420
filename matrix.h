#ifndef MATRIX_H
#define MATRIX_H
#define INDEX(A,i,j) A->cols*i+j
#define ACCESS(A,i,j) A->data[INDEX(A,i,j)]
#include<mpi.h>
#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include<string.h>
typedef struct 
{
	int rows, cols;
	double* data;	
} matrix;

void createMatrix(matrix *A)
{
	int i,j;
	A->data = malloc(A->rows*A->cols*sizeof(double));
	for(i=0;i < A->rows;i++)
	{	
		for(j=0;j < A->cols;j++)
		{
			ACCESS(A,i,j) = (double)rand() / (double)RAND_MAX;
		}
	}
}

void AllocateMatrix(matrix *A)
{
	int i,j;
	A->data = malloc(A->rows*A->cols*sizeof(double));
	for(i=0;i < A->rows;i++)
	{	
		for(j=0;j < A->cols;j++)
		{
			ACCESS(A,i,j) = 0;
		}
	}
}

void createIdentityMatrix(matrix *A)
{
	int i,j;
	A->data = malloc(A->rows*A->cols*sizeof(double));
	for(i=0; i < A->rows; i++)
	{
		for(j=0; j < A->cols; j++)
		{
			if(i == j)
			{
				ACCESS(A,i,j) = 1;
			}
			else
			{
				ACCESS(A,i,j) = 0;
			}
		}
		puts("");
	}
}

void printMatrix(matrix *A)
{
	int i,j;
	for(i=0;i < A->rows;i++)
	{	
		for(j=0;j < A->cols;j++)
		{
			printf("%0.2f ", ACCESS(A,i,j));
		}
		puts("");
	}
}

void MPIAddMatrix(matrix *A, matrix *B, matrix *C, MPI_Comm* world, int worldSize, int myRank)
{
	int N = A->rows * A->cols; //Area of amtrix
	int localLen = N/worldSize; //Length of each partition
	double *localArrA = malloc(localLen*sizeof(double));
	double *localArrB = malloc(localLen*sizeof(double));
	double *sol = malloc(localLen*sizeof(double)); // solution
	//Scatters both Matrices to be operated on in parts.
	MPI_Scatter(A->data, localLen, MPI_DOUBLE, localArrA, localLen, MPI_DOUBLE, 0, *world);
	MPI_Scatter(B->data, localLen, MPI_DOUBLE, localArrB, localLen, MPI_DOUBLE, 0, *world);
	int i;
	//Cretaes the array of solutions. The actual addition part
	//puts("Matrix A: ");
	for(i=0; i < localLen; i++)
	{
		sol[i] = localArrA[i] + localArrB[i]; 
		//printf("%0.2f ",sol[i]); 
	}
	// Puts all data into fin to be transfered to C->dataMPI_Gather(sol, localLen, MPI_DOUBLE, C->data, localLocal, MPI_DOUBLE, 0, *world);
	MPI_Gather(sol, localLen, MPI_DOUBLE, C->data, localLen, MPI_DOUBLE, 0, *world);

	free(localArrA);
	free(localArrB);
	free(sol);
}

void MPISubMatrix(matrix *A, matrix *B, matrix *C, MPI_Comm* world, int worldSize, int myRank)
{
	int N = A->rows * A->cols; //Area of amtrix
	int localLen = N/worldSize; //Length of each partition
	double *localArrA = malloc(localLen*sizeof(double));
	double *localArrB = malloc(localLen*sizeof(double));
	double *sol = malloc(localLen*sizeof(double)); // solution
	//Scatters both Matrices to be operated on in parts.
	MPI_Scatter(A->data, localLen, MPI_DOUBLE, localArrA, localLen, MPI_DOUBLE, 0, *world);
	MPI_Scatter(B->data, localLen, MPI_DOUBLE, localArrB, localLen, MPI_DOUBLE, 0, *world);
	int i;
	//Cretaes the array of solutions. The actual addition part


	for(i=0; i < localLen; i++)
	{
		sol[i] = localArrA[i] - localArrB[i];  
	}
	// Puts all data into fin to be transfered to C->data
	MPI_Gather(sol, localLen, MPI_DOUBLE, C->data, localLen, MPI_DOUBLE, 0, *world);	
	free(localArrA);
	free(localArrB);
	free(sol);
}

void MPIInnerMult(matrix *A, matrix *B, matrix *C, MPI_Comm *world, int worldSize, int myRank)
{
	int N = A->rows * B->cols; //Area of amtrix
	int localLen = N/worldSize; //Length of each partition
	double displ = malloc(myRank*sizeof(double));
	double *localArrA = malloc(localLen*sizeof(double));
	double *localArrB = malloc(localLen*sizeof(double));
	double *sol = malloc(localLen*sizeof(double)); // solution
	//Scatters both Matrices to be operated on in parts.
	//MPI_Scatterv(A->data, localLen, displ, MPI_DOUBLE, localArrA, localLen, MPI_DOUBLE, 0, *world);
	//MPI_Scatterv(B->data, localLen, displ, MPI_DOUBLE, localArrB, localLen, MPI_DOUBLE, 0, *world);
	MPI_Scatter(A->data, localLen, MPI_DOUBLE, localArrA, localLen, MPI_DOUBLE, 0, *world);
	MPI_Scatter(B->data, localLen, MPI_DOUBLE, localArrB, localLen, MPI_DOUBLE, 0, *world);
	
	int i;
	for(i=0; i < localLen; i++)
	{
		sol[i] += A->data[i] * B->data[i]; 
		//printf("Sol at %d: %d\n", i, sol[i]);
	}
	//MPI_Gatherv(sol, localLen, MPI_DOUBLE, C->data, localLen, displ, MPI_DOUBLE, 0, *world);	
	MPI_Gather(sol, localLen, MPI_DOUBLE, C->data, localLen, MPI_DOUBLE, 0, *world);	

	MPI_FILE fh;
	int offset = N*myRank*sizeof();
	char buf[N];
	int i,j;
	for(i=0;i < A.rows; i++)
	{
		for(j=0;j<A.cols;j++)
		{
			buf[i] = ACCESS(A,i,j);
		}
	}
	MPI_File_open(world,"matrix.txt",MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
	MPI_File_write(fh, buf, 0, MPI_DOUBLE, MPI_STATUS_IGNORE);
	MPI_File_close(&fh);

	free(localArrA);
	free(localArrB);
	free(sol);
}
/*
void BlockedMatrixMult(matrix *A, matrix *B, matrix *C, MPI_Comm *world, int worldSize, int myRank, int bSize)
{
	int i,j,k,l,m;
	double sum;
	int n = A->rows;
	int en = bSize * (n/bSize);
	for(i=0;i<en;i+=bSize)
	{
		for(j=0;j<en;j+=bSize)
		{
			for(k=0;k<n;k++)
			{
				for(l=j;l<j+bSize;l++)
				{
					sum = ACCESS(C,i,j);
					for(m=i;m<i+bSize;m++)
					{
						sum+= A->data[m] * B->data[m];	
					}
					C->data = sum;
				}
			}
		}
	}			
}
*/
//Row reduction formula
//Can compute inverse of a matrix
void MPIGaussJordan(matrix *A, matrix *B, matrix *D, MPI_Comm *world, int worldSize, int myRank)
{
	//puts("Hi");
	int n = A->rows * A->cols; //Area of amtrix
	int i,j,k;
	char L[n];	
	int localLen = n/worldSize; //Length of each partition
	double *localArrA = malloc(localLen*sizeof(double));
	double *localArrB = malloc(localLen*sizeof(double));
	double *sol = malloc(localLen*sizeof(double));
	//Scatters both Matrices to be operated on in parts.

	for(i=0; i < n; i++)
	{
		L[i] = i % worldSize;
	}
	
	MPI_Scatter(A->data, localLen, MPI_DOUBLE, localArrA, localLen, MPI_DOUBLE, 0, *world);
	MPI_Scatter(B->data, localLen, MPI_DOUBLE, localArrB, localLen, MPI_DOUBLE, 0, *world);
	//K is the pivot row
	for(k=0;k<localLen;k++)
	{
		for(i=k+1; i<localLen; i++)
		{	
			if(L[k] == myRank)
			{
				sol[i] = localArrA[i] / localArrB[i];
			}
		}

		for(i=k+1; i<localLen; i++)
		{
			if(L[i] == myRank)
			{
				ACCESS(A,i,j) = ACCESS(A,i,j) - (sol[i] * ACCESS(A,j,k));
			}
		 	ACCESS(B,i,j) = ACCESS(D,i,j) - (sol[i] * ACCESS(B,j,k));
		}	
	}
	MPI_Gather(sol, localLen, MPI_DOUBLE, D->data, localLen, MPI_DOUBLE, 0, *world);	
	free(localArrA);
	free(localArrB);
	free(sol);
}

#endif
