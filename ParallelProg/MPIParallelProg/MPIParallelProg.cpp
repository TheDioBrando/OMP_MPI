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
	int sumCount[2];
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

	sumCount[0] = 0, sumCount[1] = 0;
	for (int i = 0; i < partition; i++) {
		if (aPart[i] > 0) {
			sumCount[1]++;
			sumCount[0] += aPart[i];
		}
	}

	MPI_Reduce(sumCount, sumCount, 2, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

	if (!rank)
	{
		double avg = (double)sumCount[0] / sumCount[1];
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
	int rank, size, localmaxmin, localminmax, maxmin, minmax, localmin, localmax;
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

	for (int i = 0; i < partition; i++)
	{
		localmin = part[columns * i];
		localmax = part[columns * i];
		for (int j = 0; j < columns; j++)
		{
			if (localmax < part[columns * i + j])
			{
				localmax = part[columns * i + j];
			}

			if (localmin > part[columns * i + j])
			{
				localmin = part[columns * i + j];
			}
		}
		if (i == 0)
		{
			localminmax = localmax;
			localmaxmin = localmin;
		}
		if (localmaxmin < localmin)
		{
			localmaxmin = localmin;
		}
		if (localminmax > localmax)
		{
			localminmax = localmax;
		}
	}

	MPI_Reduce(&localmaxmin, &maxmin, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
	MPI_Reduce(&localminmax, &minmax, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
	if (!rank)
	{
		printf("Maxmin = %d, Minmax = %d\n", maxmin, minmax);
		if (maxmin == minmax)
		{
			printf("Maxmin = Minmax\n");
		}
		else
		{
			printf("Maxmin != Minmax\n");
		}
	}

	MPI_Finalize();
}

void task8()
{
	const int N = 10;
	int rank, size;
	int a[N];
	int recvA[N];

	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int partition = (int)round((double)N / size);
	int* part = new int[partition];

	if (rank == 0)
	{
		for (int i = 0; i < N; i++)
		{
			a[i] = rand() % 100;
			printf("a root = %d\n", a[i]);
		}

		for (int i = 0; i < size; ++i) {
			int* sendPart = new int[partition];
			for (int j = 0; j < partition; j++)
			{
				sendPart[j] = a[i * partition + j];
			}
			if (i == 0)
			{
				part = sendPart;
			}
			else
			{
				MPI_Send(sendPart, partition, MPI_INT, i, i, MPI_COMM_WORLD);
			}
		}
	}

	if (rank)
	{
		MPI_Recv(part, partition, MPI_INT, 0, rank, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
	}
	for (int i = 0; i < partition; i++)
	{
		printf("rank = %d, %d\n",rank, part[i]);
	}

	if (rank) {
		MPI_Send(part, partition, MPI_INT, 0, rank, MPI_COMM_WORLD);
	}
	else
	{
		for (int i = 0; i < partition; i++) {
			recvA[i] = part[i];
		}
	}

	if (!rank)
	{
		for (int i = 1; i < size; i++)
		{
			MPI_Recv(part, partition, MPI_INT, i, i, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
			for (int j = 0; j < partition; j++)
			{
				recvA[i * partition + j] = part[j];
			}
		}
		for (int i = 0; i < N; i++)
			printf("result arr = %d\n", recvA[i]);
	}

	MPI_Finalize();
}

void task11()
{
	int rank, size;
	int msg = 3;
	int send, recv;
	int sendVal, recvVal;

	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	send = rank + 1;
	recv = rank - 1;
	if (!rank)
		recv = size - 1;
	if (rank == size - 1)
		send = 0;

	if (!rank)
	{
		sendVal = msg * 3;
		MPI_Send(&sendVal, 1, MPI_INT, send, send, MPI_COMM_WORLD);
		printf("proc %d, send msg %d, to %d\n", rank, sendVal, send);
	}
	else
	{
		MPI_Recv(&recvVal, 1, MPI_INT, recv, rank, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
		printf("proc %d, receive msg %d from %d\n", rank, recvVal, recv);
		sendVal = recvVal * 3;
		MPI_Send(&sendVal, 1, MPI_INT, send, send, MPI_COMM_WORLD);
		printf("proc %d, send msg %d, to %d\n", rank, sendVal, send);
	}

	if (!rank)
	{
		MPI_Recv(&recvVal, 1, MPI_INT, recv, rank, MPI_COMM_WORLD, MPI_STATUSES_IGNORE);
		printf("proc %d, receive msg %d from %d\n", rank, recvVal, recv);
	}

	MPI_Finalize();
}

void generateArray(int* arr, int size)
{
	srand(time(NULL));

	for (int i = 0; i < size; i++)
	{
		arr[i] = rand() % 100;
	}
}

void LocalQuickSort(int* arr, int start, int end)
{
	if (start == end)
		return;

	int pivot = arr[end];
	int storedIndex = start;

	for (int i = start; i < end; i++)
	{
		if (arr[i] <= pivot)
		{
			int temp = arr[storedIndex];
			arr[storedIndex] = arr[i];
			arr[i] = temp;
			storedIndex = storedIndex + 1;
		}
	}

	int temp = arr[storedIndex];
	arr[storedIndex] = arr[end];
	arr[end] = temp;

	if (storedIndex > start)
		LocalQuickSort(arr, start, storedIndex - 1);
	if (storedIndex < end)
		LocalQuickSort(arr, storedIndex + 1, end);
}

void PivotDistribution(int* part, int partSize, int Dim,
	int Mask, int Iter, int* pPivot) {
	int rank, ProcNum;
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	MPI_Group WorldGroup;
	MPI_Group SubcubeGroup; // группа процессов - подгиперкуб
	MPI_Comm SubcubeComm; // коммуникатор подгиперкуба
	int j = 0;
	int GroupNum = ProcNum / (int)pow(2, Dim - Iter);
	int* ProcRanks = new int[GroupNum];
	// формирование списка рангов процессов для гиперкуба
	int StartProc = rank - GroupNum;
	if (StartProc < 0) StartProc = 0;
	int EndProc = rank + GroupNum;
	if (EndProc > ProcNum) EndProc = ProcNum;
	for (int proc = StartProc; proc < EndProc; proc++) {
		if ((rank & Mask) >> (Iter) == (proc & Mask) >> (Iter)) {
			ProcRanks[j++] = proc;
		}
	}
	//объединение процессов подгиперкуба в одну группу
	MPI_Comm_group(MPI_COMM_WORLD, &WorldGroup);
	MPI_Group_incl(WorldGroup, GroupNum, ProcRanks, &SubcubeGroup);
	MPI_Comm_create(MPI_COMM_WORLD, SubcubeGroup, &SubcubeComm);
	// поиск и рассылка ведущего элемента всем процессам подгиперкуба
	if (rank == ProcRanks[0])
		*pPivot = part[(partSize) / 2];
	MPI_Bcast(pPivot, 1, MPI_DOUBLE, 0, SubcubeComm);
	MPI_Group_free(&SubcubeGroup);
	MPI_Comm_free(&SubcubeComm);
	delete[] ProcRanks;
}

int GetProcDataDivisionPos(int* part, int size, int pivot)
{
	if (part[0] > pivot)
		return 0;
	if (part[size - 1] < pivot)
		return size;

	int left = 0, right = size - 1;

	while (left < right)
	{
		int middle = (right + left) / 2;
		if (pivot <= part[middle])
			right = middle;
		else
			left = middle + 1;
	}
	
	return right;
}

void DataMerge(int* pMergeData, int MergeDataSize, int* pData, int DataSize,
	int* pRecvData, int RecvDataSize)
{
	for (int i = 0; i < DataSize; i++)
		pMergeData[i] = pData[i];
	
	for (int i = 0; i < RecvDataSize; i++)
		pMergeData[DataSize + i] = pRecvData[i];
}

void ParallelQuickSort(int* part, int size)
{
	int rank, ProcNum;
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);


	MPI_Status status;
	int CommProcRank; // ранг процессора, с которым выполняется взаимодействие
	int* pData, // часть блока, остающаяся на процессоре
		* pSendData, // часть блока, передаваемая процессору CommProcRank
		* pRecvData, // часть блока, получаемая от процессора CommProcRank
		* pMergeData; // блок данных, получаемый после слияния
	int DataSize, SendDataSize, RecvDataSize, MergeDataSize;
	int HypercubeDim = (int)ceil(log(ProcNum) / log(2)); // размерность гиперкуба
	int Mask = ProcNum;
	int Pivot;

	LocalQuickSort(part, 0, size - 1);
	// итерации обобщенной быстрой сортировки
	for (int i = HypercubeDim; i > 0; i--) {
		// определение ведущего значения и его рассылка всем процессорам
		PivotDistribution(part, size, HypercubeDim, Mask, i, &Pivot);
		Mask = Mask >> 1;
		// определение границы разделения блока
		int pos = GetProcDataDivisionPos(part, size, Pivot);
		// разделение блока на части
		if (((rank & Mask) >> (i - 1)) == 0) { // старший бит = 0
			pSendData = &part[pos + 1];
			SendDataSize = size - pos- 1;
			if (SendDataSize < 0) 
				SendDataSize = 0;
			CommProcRank = rank + Mask;
			pData = &part[0];
			DataSize = pos + 1;
		}
		else { // старший бит = 1
			pSendData = &part[0];
			SendDataSize = pos + 1;
			if (SendDataSize > size)
				SendDataSize = pos;
			CommProcRank = rank - Mask;
			pData = &part[pos + 1];
			DataSize = size - pos - 1;
			if (DataSize < 0) DataSize = 0;
		}
		// пересылка размеров частей блоков данных
		MPI_Sendrecv(&SendDataSize, 1, MPI_INT, CommProcRank, 0,
			&RecvDataSize, 1, MPI_INT, rank, 0, MPI_COMM_WORLD, &status);
		// пересылка частей блоков данных
		pRecvData = new int[RecvDataSize];
		MPI_Sendrecv(pSendData, SendDataSize, MPI_DOUBLE,
			CommProcRank, 0, pRecvData, RecvDataSize, MPI_DOUBLE,
			rank, 0, MPI_COMM_WORLD, &status);
		// слияние частей
		MergeDataSize = DataSize + RecvDataSize;
		pMergeData = new int[MergeDataSize];
		DataMerge(pMergeData, MergeDataSize, pData, DataSize, pRecvData, RecvDataSize);
		delete[] part;
		delete[] pRecvData;
		part = pMergeData;
		size = MergeDataSize;
	}
}

void Sort()
{
	int N = 100;
	int rank, size;
	int* arr = new int[N];

	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int partition = N / size;
	int* part = new int[partition];

	int* ind = new int[size];
	int* len = new int[size];

	printf("Hello");

	if (!rank)
	{
		generateArray(arr, N);
		//LocalQuickSort(arr, 0, N - 1);

		for (int i = 0; i < N; i++)
		{
			printf("arr = %d\n", arr[i]);
		}

		for (int i = 0; i < size; ++i) {
			ind[i] = i * partition;
			len[i] = partition;
		}
	}

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Scatterv(arr, len, ind, MPI_INT, part, partition, MPI_INT, 0, MPI_COMM_WORLD);

	ParallelQuickSort(part, partition);
	partition = sizeof(part);
	MPI_Gather(&partition, 1, MPI_INT, len, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if (!rank)
	{
		ind[0] = 0;
		for (int i = 1; i < size; i++)
		{
			ind[i] = ind[i - 1] + len[i - 1];
		}
	}

	MPI_Gatherv(part, partition, MPI_INT, arr, len, ind, MPI_INT, 0, MPI_COMM_WORLD);

	if (!rank)
	{
		printf("sorted array \n");
		for (int i = 0; i < N; i++)
		{
			printf("arr = %d\n", arr[i]);
		}
	}
	delete[] arr;
	delete[] part;
	delete[] ind;
	delete[] len;

	MPI_Finalize();
}

int main()
{
	printf("Heloo");
	Sort();
}
