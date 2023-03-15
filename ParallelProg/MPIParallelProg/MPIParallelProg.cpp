#include <iostream>
#include <mpi.h>

void task1()
{
    int rank, size;

    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    printf("Hello from rank %d, size %d", rank, size);

    MPI_Finalize();
}

void task2()
{
    const int N = 100;
    int rank, size;
    int a[N];

    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (rank == 0) 
	{
		for (int i = 0; i < N; i++)
		{
			a[i] = rand() % 100;
		}
	}

	MPI_Bcast(a, N, MPI_INT, 0, MPI_COMM_WORLD);

	int k = N / size;
	int i1 = k * rank;
	int i2 = k * (rank + 1);
	if (rank == size - 1)
		i2 = N;

	int max = a[i1], totalMax = a[0];
	for (int i = i1 + 1; i < i2; i++)
		if (max < a[i])
			max = a[i];

	printf("max = %d, rank = %d\n", max, rank);
	MPI_Reduce(&max, &totalMax, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);

	if (!rank)
	{
		printf("totalMax = %d\n", totalMax);
	}

    MPI_Finalize();
}

void task3()
{
	int rank, size;
	int count = 1000000;
	double x, y;

	int dots_count = 0, local_dots_count = 0;
	MPI_Init(NULL, NULL);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	srand(time(NULL));
	int partition = count / size;

	for (int i = 0; i < partition; ++i) {
		x = (double)rand() / (double)RAND_MAX * 2 - 1;
		y = (double)rand() / (double)RAND_MAX * 2 - 1;

		if (x*x + y*y <= 1)
			local_dots_count++;
	}
	printf("%d \n", local_dots_count);
	MPI_Reduce(&local_dots_count, &dots_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

	if (rank == 0) {
		double answer = (double)(4 * dots_count) / (double)count;
		printf("Pi = %f    %d", answer, dots_count);
	}
	MPI_Finalize();
}

void task4()
{
	const int N = 100;
	int rank, size, totalCount = 0, count = 0, totalSum = 0, sum = 0;
	int a[N];

	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int partition = (int)round((double)N / size);
	if (N % size != 0) {
		if (rank == size - 1) partition = N - partition * (size - 1);
	}

	int* aPart = new int[partition];

	int* ind = new int[size];
	int* len = new int[size];


	if (rank == 0)
	{
		for (int i = 0; i < N; i++)
		{
			if (rand() % 2 == 0)
			{
				a[i] = rand() % 100;
			}
			else
			{
				a[i] = -1* rand() % 100;
			}
		}
		
		for (int i = 0; i < size; ++i) {
			ind[i] = i * partition;
			len[i] = partition;
		}

		if (N % size != 0) {
			len[size - 1] = N - partition * (size - 1);
		}
	}
	
	MPI_Scatterv(a, len, ind, MPI_INT, aPart, partition, MPI_INT, 0, MPI_COMM_WORLD);

	for (int i = 0; i < partition; i++) {
		if (aPart[i] > 0) {
			count++;
			sum += aPart[i];
		}
	}

	MPI_Reduce(&count, &totalCount, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&sum, &totalSum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

	if (!rank)
	{
		double avg = (double)totalSum / totalCount;
		printf("Avg = %f", avg);
	}

	MPI_Finalize();
}

void task5() 
{
	const int N = 100;
	int rank, size, sum = 0, totalSum = 0;
	int a[N], b[N];

	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int partition = (int)round((double)N / size);
	if (N % size != 0) {
		if (rank == size - 1) partition = N - partition * (size - 1);
	}


	int* aPart = new int[partition];
	int* bPart = new int[partition];

	int* ind = new int[size];
	int* len = new int[size];

	if (rank == 0)
	{
		for (int i = 0; i < N; i++)
		{
			a[i] = rand() % 100;
			b[i] = rand() % 100;
		}

		for (int i = 0; i < size; ++i) {
			ind[i] = i * partition;
			len[i] = partition;
		}

		if (N % size != 0) {
			len[size - 1] = N - partition * (size - 1);
		}
	}

	MPI_Scatterv(a, len, ind, MPI_INT, aPart, partition, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Scatterv(b, len, ind, MPI_INT, bPart, partition, MPI_INT, 0, MPI_COMM_WORLD);

	for (int i = 0; i < partition; i++)
	{
		sum += aPart[i] * bPart[i];
	}

	MPI_Reduce(&sum, &totalSum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	if (!rank)
		printf("Scalar of a and b = %d", totalSum);

	MPI_Finalize();
}

void task6()
{
	const int N = 100;
	int columns = N, rows = N;
	int rank, size, min, max, totalMax, totalMin;
	int** matrix = new int*[N];
	for (int i = 0; i < N; i++)
	{
		matrix[i] = new int[N];
	}

	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int partition = (int)round((double)N / size);
	if (N % size != 0) {
		if (rank == size - 1) partition = N - partition * (size - 1);
	}

	int* ind = new int[size];
	int* len = new int[size];
	int* sendArr = new int[N * N];
	int* part = new int[partition * columns];

	if (rank == 0)
	{
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++)
			{
				matrix[i][j] = rand() % 100;
				sendArr[columns * i + j] = matrix[i][j];
			}
		}

		for (int i = 0; i < size; ++i) {
			ind[i] = i * partition*columns;
			len[i] = partition*columns;
		}

		if (rows % size != 0) {
			len[size - 1] = (rows - partition * (size - 1))*columns;
		}
	}

	MPI_Scatterv(sendArr, len, ind, MPI_INT, part, partition * columns, MPI_INT, 0, MPI_COMM_WORLD);

	max = part[0], min = part[0];

	for (int i = 0; i < partition; i++)
	{
		for (int j = 0; j < columns; j++)
		{
			if (max < part[columns * i + j])
			{
				max = part[columns * i + j];
			}

			if (min > part[columns * i + j])
			{
				min = part[columns * i + j];
			}
		}
	}

	MPI_Reduce(&min, &totalMin, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
	MPI_Reduce(&max, &totalMax, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
	if (!rank)
	{
		printf("Max = %d, Min = %d", totalMax, totalMin);
	}

	MPI_Finalize();
}

int main()
{
    task6();
}
