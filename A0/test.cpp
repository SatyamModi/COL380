#include <bits/stdc++.h>
#include <omp.h>
#include <mpi.h>
using namespace std;

int main()
{
	// int A[4][4];
	// int B[4][4];
	// int rank, size;

	// if (rank == 0)
	// {
	// 	for (int i = 0; i < 4; i++)
	// 	{
	// 		for (int j = 0; j < 4; j++)
	// 		{
	// 			A[i][j] = 8;
	// 		}
	// 	}
	// }
	// else
	// {
	// 	for (int i = 0; i < 4; i++)
	// 	{
	// 		for (int j = 0; j < 4; j++)
	// 		{
	// 			A[i][j] = 5;
	// 		}
	// 	}
	// }

	// MPI_Init(NULL, NULL);
	// MPI_Comm_rank(MPI_COMM_WORLD, &rank);
 //    MPI_Comm_size(MPI_COMM_WORLD, &size);

 //    int recvcounts[1] = {8};
 //    int displs[1] = {8};
 //    if (rank == 0)
	// 	MPI_Gatherv(&(A[0][0]), 8, MPI_INT, &(B[0][0]), recvcounts, displs, MPI_INT, 0, MPI_COMM_WORLD);
	// else
	// 	MPI_Gatherv(&(A[0][0]), 8, MPI_INT, NULL, NULL, NULL, MPI_INT, 0, MPI_COMM_WORLD);
	// if (rank == 0)
	// {
	// 	for (int i = 0; i < 4; i++)
	// 	{
	// 		for (int j = 0; j < 4; j++)
	// 		{
	// 			cout << B[i][j] << "\n";
	// 		}
	// 	}
	// }
	// MPI_Finalize();
	MPI_Init(NULL, NULL);

	int size, rank;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int root = 0;
	int B[16];
	if (rank == 0)
	{
		int my_val[1][4] = {{4, 4, 4, 4}};
		int counts[3] = {4,4,8};
		int displs[3] = {0,4,8};
		
		MPI_Gatherv(&my_val, 4, MPI_INT, B, counts, displs, MPI_INT, root, MPI_COMM_WORLD);
	}
	else if (rank == 1)
	{
		int my_val[1][4] = {{5, 5, 5, 5}};
		MPI_Gatherv(&(my_val[0][0]), 4, MPI_INT, NULL, NULL, NULL, MPI_INT, root, MPI_COMM_WORLD);
	}

	else
	{
		int my_val[2][4] = {{6, 6, 6, 6}, {7,7,7,7}};
		MPI_Gatherv(&(my_val[0][0]), 8, MPI_INT, NULL, NULL, NULL, MPI_INT, root, MPI_COMM_WORLD);
	}

	if (rank == 0)
	{
		for (int i = 0; i < 16; i++)
		{
			cout << B[i] << endl;
		}
	}
	MPI_Finalize();
	return 0;
}