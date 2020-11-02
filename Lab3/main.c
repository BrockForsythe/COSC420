#include "matrix.h"

int main(int argc, char **argv)
{
	//Initialize MPI
	MPI_Init(&argc,&argv);
	MPI_Comm world = MPI_COMM_WORLD;
	int myRank, worldSize;
	MPI_Comm_rank(world,&myRank);	
	MPI_Comm_size(world,&worldSize);
	//Allows for random numbers
	srand(time(NULL) + myRank);
	//Initialze the matrices from struct matrix
	matrix A;
	matrix B; 
	matrix C;
	matrix D;
	matrix E;
	A.rows = atoi(argv[1]);
	A.cols = atoi(argv[2]);
	B.rows = atoi(argv[3]);
	B.cols = atoi(argv[4]);
	C.rows = A.rows;
	C.cols = B.cols;
	D.rows = C.rows;
	D.cols = C.cols;
	E.rows = A.rows;
	E.cols = A.cols;
	//Fills A and B with random nums and C with an empty arra

	createMatrix(&A);
	createMatrix(&B);
	AllocateMatrix(&C);
	AllocateMatrix(&D);	
	createIdentityMatrix(&E);	
	//Number of nodes to use in the matrix multiplication
	int N = worldSize;
	//int bSize = N;	
	//Keeps N a whole number
	if(N % worldSize > 0)
	{
		printf("Choose N divisible by %d\n", worldSize);
		MPI_Finalize();
		return 0;
	}
	
	//Prints matrices
	if(myRank == 0)
	{	
		printf("\nPrinting Matrix A \n");
		printMatrix(&A);
	}
	if(myRank == 0)	
	{
		printf("\nPrinting Matrix B \n");
		printMatrix(&B);
	}
	//int bSize = N/worldSize;
	//BlockedMatrixMult(&A,&B,&C,&world,worldSize,myRank,bSize);
	MPIInnerMult(&A,&B,&C, &world, worldSize, myRank);
	if(myRank == 0)
	{
		printf("\nPrinting multiplied matrix C \n");
		printMatrix(&C);
	}
	
	if(myRank == 0)
	{
		printf("\nPrinting Identity Matrix \n");
		printMatrix(&E);
	}

	if((A.rows == B.rows) && (A.cols == B.cols))
	{	
		MPIGaussJordan(&A, &E, &D, &world, worldSize, myRank);
	}
	
	if(myRank == 0)
	{
		printf("\nPrinting inverted matrix \n");
		printMatrix(&D);
	}	

	MPI_Finalize();
	return 0; 
}
