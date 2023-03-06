// ParallelProg.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <ctime>
#include <omp.h>

using namespace std;

void task1()
{
    #pragma omp parallel num_threads(8)
        {
            int rank = omp_get_thread_num();
            int size = omp_get_num_threads();
            printf("Thread %d out of %d threads. Hello world!\n", rank, size);
        }
}

void task2()
{
    int thread_num = 3;
    omp_set_num_threads(thread_num);
    
#pragma omp parallel if(thread_num > 2) 
    {
        int rank = omp_get_thread_num();
        int size = omp_get_num_threads();
        printf("Thread %d out of %d threads parallel 1\n", rank, size);
    }

    thread_num = 2;
    omp_set_num_threads(thread_num);
#pragma omp parallel if(thread_num > 2) 
    {
        int rank = omp_get_thread_num();
        int size = omp_get_num_threads();
        printf("Thread %d out of %d threads parallel 2\n", rank, size);
    }
}

void task3()
{
    int a = 1, b = 3;

    printf("Sequential before parallel a = %d, b = %d \n",a ,b);

#pragma omp parallel num_threads(2) private(a) firstprivate(b)
    {
        int rank = omp_get_thread_num();
        a += rank; b += rank;

        printf("Parallel a = %d, b = %d \n", a, b);
    }

    printf("Sequential after parallel a = %d, b = %d\n", a, b);

    printf("Sequential before parallel a = %d, b = %d\n", a, b);

#pragma omp parallel num_threads(4) shared(a) private(b)
    {
        int rank = omp_get_thread_num();
        a -= rank; b -= rank;

        printf("Parallel a = %d, b = %d\n", a, b);
    }

    printf("Sequential after parallel a = %d, b = %d\n", a, b);
}

void task4()
{
    int* a = new int[10];
    int* b = new int[10];

    for (int i = 0; i < 10; i++)
    {
        a[i] = i; b[i] = i * 2;
    }

#pragma omp parallel num_threads(2)
    {
        int rank = omp_get_thread_num();

        if (rank == 0)
        {
            int min = a[0];
            for (int i = 1; i < 10; i++)
            {
                if (a[i] < min)
                    min = a[i];
            }

            printf("min a = %d \n", min);
        }
        else
        {
            int max = b[0];
            for (int i = 1; i < 10; i++)
            {
                if (b[i] > max)
                    max = b[i];
            }

            printf("max b = %d \n", max);
        }
    }
}

void task5()
{
    srand(time(NULL));

    int** d = new int* [6];
    for (int i = 0; i < 6; i++)
    {
        d[i] = new int[8];
    }

    for (int i = 0; i < 6; i++)
    {
        for (int j = 0; j < 8; j++)
        {
            d[i][j] = rand() % 100;
        }
    }

#pragma omp parallel
    {
#pragma omp sections
        {
#pragma omp section
            {
                int sum = 0;
                int count = 6*8;
                for (int i = 0; i < 6; i++)
                {
                    for (int j = 0; j < 8; j++)
                    {
                        sum += d[i][j];
                    }
                }

                double result = sum * 1.0 / count;
                printf("Section - 0 thread - %d, mean = %f\n", omp_get_thread_num(), result);
            }
#pragma omp section
            {
                int max = d[0][0];
                int min = d[0][0];
                for (int i = 0; i < 6; i++)
                {
                    for (int j = 0; j < 8; j++)
                    {
                        if (max < d[i][j])
                            max = d[i][j];
                        if (min > d[i][j])
                            min = d[i][j];
                    }
                }
                printf("Section - 1 thread - %d, min = %d, max = %d\n", omp_get_thread_num(), min, max);
            }
#pragma omp section
            {
                int count = 0;
                for (int i = 0; i < 6; i++)
                {
                    for (int j = 0; j < 8; j++)
                    {
                        if (d[i][j] % 3 == 0)
                            count++;
                    }
                }
                printf("Section - 2 thread - %d, count = %d\n", omp_get_thread_num(), count);
            }
        }
    }

}

void task6()
{
    int* a = new int[100];
    for (int i = 0; i < 100; i++)
        a[i] = i;

    int sum = 0;
#pragma omp parallel for reduction(+ : sum)
    for (int i = 0; i < 100; i++)
    {
        sum += a[i];
    }

    double mean = sum * 1.0 / 100;
    printf("reduction mean = %f \n", mean);

    sum = 0;

#pragma omp parallel for shared(sum)
    for (int i = 0; i < 100; i++)
    {
        sum += a[i];
    }

    mean = sum * 1.0 / 100;
    printf("without reduction mean = %f \n", mean);
}

void task7()
{
    int a[12], b[12], c[12];
#pragma omp parallel for num_threads(3) schedule(static, 4)
    for (int i = 0; i< 12; i++)
    {
        a[i] = i; b[i] = i * 2;
        int rank = omp_get_thread_num();
        int size = omp_get_num_threads();
        printf("rank = %d, size = %d, a[%d] = %d, b[%d]=%d\n", rank, size, i, a[i], i, b[i]);
    }

#pragma omp parallel for num_threads(4) schedule(dynamic, 2)
    for (int i = 0;i<12;i++)
    {
        c[i] = a[i] + b[i];
        int rank = omp_get_thread_num();
        int size = omp_get_num_threads();
        printf("rank = %d, size = %d, c[%d] = %d\n", rank, size, i, c[i]);
    }
}

void task8()
{
    int arr[16000];
    double b[16000];
    int chunk = 40;
#pragma omp parallel for schedule(dynamic, 4)
    for (int i = 0; i < 16000; i++)
    {
        arr[i] = i;
    }

    double time = omp_get_wtime();

#pragma omp parallel for num_threads(8) schedule(static, chunk)
    for (int i = 1; i < 15999; i++)
    {
        b[i] = (arr[i - 1] + arr[i] + arr[i + 1]) / 3.0;
    }
    printf("time to static = %f\n", omp_get_wtime() - time);

    time = omp_get_wtime();

#pragma omp parallel for num_threads(8) schedule(dynamic, chunk)
    for (int i = 1; i < 15999; i++)
    {
        b[i] = (arr[i - 1] + arr[i] + arr[i + 1]) / 3.0;
    }
    printf("time to dynamic = %f\n", omp_get_wtime() - time);

    time = omp_get_wtime();

#pragma omp parallel for num_threads(8) schedule(guided, chunk)
    for (int i = 1; i < 15999; i++)
    {
        b[i] = (arr[i - 1] + arr[i] + arr[i + 1]) / 3.0;
    }
    printf("time to guided = %f\n", omp_get_wtime() - time);

    time = omp_get_wtime();

#pragma omp parallel for num_threads(8) schedule(runtime)
    for (int i = 1; i < 15999; i++)
    {
        b[i] = (arr[i - 1] + arr[i] + arr[i + 1]) / 3.0;
    }
    printf("time to runtime = %f\n", omp_get_wtime() - time);
}

void task9()
{
    int n = 4000;
    int** matrix = new int* [n];
    int* vec = new int[n];
    int* result = new int[n];
    for (int i = 0; i < n; i++)
        matrix[i] = new int[n];

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
            matrix[i][j] = i + j;
        vec[i] = i;
    }

    double time = omp_get_wtime();

    int sum;
    for (int i = 0; i < n; i++)
    {
        sum = 0;
        for (int j = 0; j < n; j++)
        {
            sum += matrix[i][j] * vec[j];
        }
        result[i] = sum;
    }

    printf("seq time = %f, result = %d\n", omp_get_wtime() - time, result[n-2]);

    result = new int[n];
    time = omp_get_wtime();

#pragma omp parallel for schedule(guided, 400) private(sum)
    for (int i = 0; i < n; i++)
    {
        sum = 0;
        for (int j = 0; j < n; j++)
        {
            sum += matrix[i][j] * vec[j];
        }
        result[i] = sum;
    }

    printf("parallel time = %f, result = %d\n", omp_get_wtime() - time, result[n - 2]);
}

void task10()
{
    srand(time(NULL));

    int** d = new int* [6];
    for (int i = 0; i < 6; i++)
    {
        d[i] = new int[8];
    }

    for (int i = 0; i < 6; i++)
    {
        for (int j = 0; j < 8; j++)
        {
            d[i][j] = rand() % 100;
        }
    }

    int min = d[0][0], max = d[0][0];
#pragma omp parallel for schedule(dynamic, 2) num_threads(4)
    for (int i = 0; i < 6; i++)
    {
        for (int j = 0; j < 8; j++)
        {
            if (min > d[i][j])
            {
#pragma omp critical
                {
                    if (min > d[i][j])
                        min = d[i][j];
                }
            }
            if (max < d[i][j])
            {
#pragma omp critical
                {
                    if (max < d[i][j])
                        max = d[i][j];
                }
            }
        }
    }

    printf("max = %d, min = %d", max, min);
}

void task11()
{
    srand(time(NULL));

    int* a = new int [30];

    for (int i = 0; i < 30; i++)
    {
        a[i] = rand() % 100;
    }

    int count = 0;
#pragma omp parallel for num_threads(4)
    for (int i = 0; i < 30; i++)
    {
        if (a[i] % 9 == 0)
        {
#pragma omp atomic
            count++;
        }
    }

    printf("Count = %d", count);
}

void task12()
{
    srand(time(NULL));

    omp_lock_t lock;
    int n = 50;
    int* a = new int[n];
    int max;

    for (int i = 0; i < n; i++)
    {
        a[i] = rand() % 100;
    }

    max = a[0];
    omp_init_lock(&lock);
#pragma omp parallel for
    for (int i = 1; i < n; i++)
    {
        if ((a[i]%7 == 0) 
            && max < a[i])
        {
            omp_set_lock(&lock);
            if (max < a[i])
                max = a[i];
            omp_unset_lock(&lock);
        }
    }

    printf("max = %d", max);
    omp_destroy_lock(&lock);
}

int main()
{
    //task1();
    //task2();
    //task3();
    //task4();
    //task5();
    //task6();
    //task7();
    //task8();
    //task9();
    //task10(); //to do
    //task11();
    //task12();

}
