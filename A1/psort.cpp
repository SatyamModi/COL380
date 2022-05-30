#include "psort.h"
#include <omp.h>

uint32_t Threshold;
uint32_t Threshold_set=0;

struct Partition
{
    uint32_t size;
    uint32_t *array;

    Partition(uint32_t n)
    {
        size = n;
        array = new uint32_t[size];
    }
};

void swap(uint32_t *a, uint32_t *b)
{
    uint32_t temp = *a;
    *a = *b;
    *b = temp;
}

int Partitioning(uint32_t *data, int low, int high)
{
    uint32_t pivot = data[high];
    int i = low-1;

    for (int j=low; j <= high-1; j++)
    {
        if (data[j] <= pivot)
        {
            i++;
            swap(&data[i], &data[j]);
        }
    }

    swap(&data[i+1], &data[high]);
    return (i+1);
}

void QuickSort(uint32_t *data, int low, int high)
{
    if (low < high)
    {
        int pi = Partitioning(data, low, high);
        QuickSort(data, low, pi-1);
        QuickSort(data, pi+1, high);
    }
}

void SequentialSort(uint32_t *data, uint32_t n)
{
    QuickSort(data, 0, n-1);
}

void ParallelSort(uint32_t *data, uint32_t n, int p)
{
    // Entry point to your sorting implementation.
    // Sorted array should be present at location pointed to by data.

    // If n <= Threshold
    if (!Threshold_set)
    {
        Threshold = 2*n/p;
        Threshold_set=1;
    }

    if (n <= Threshold) 
    {
        SequentialSort(data, n);
        return;
    }

    // Split the array A into A_0, A_1, ..., A_p-1 of size n/p +- 1. 
    // All A_i are contiguous
    uint32_t partition_size = n/p;
    
    int i,j;

    // selecting pseudo splitters into R where R = [r_0, r_1, .., r_p*(p-1)]
    uint32_t *R = new uint32_t[p*p];
    
    uint32_t rsize = 0;
    for (i = 0; i < p; i++)
    {
        for (j = 0; j < p; j++)
        {
            R[rsize] = data[partition_size*i+j];
            rsize += 1;
        }
    }

    SequentialSort(R, rsize);

    j = 0;
    uint32_t S[p-1];
    for (j = 0; j < p-1; j++)
    {
        S[j] = R[(j+1)*p];
    }

    // size[p] array stores the size of buckets and idx maps number A[i] to bucket B[j]
    uint32_t *BucketSize = new uint32_t[p];
    for (i = 0; i < p; i++)
    {
        BucketSize[i] = 0;
    }

    uint32_t *idx = new uint32_t[n];
    uint32_t k = 0;

    for (k = 0; k < n; k++)
    {
        idx[k] = 0;
    }

    k = 0;
    for (k = 0; k < n; k++)
    {
        if (data[k] <= S[0])
        {
            BucketSize[0]++;
            idx[k] = 0;
        }
        else if (data[k] > S[p-2])
        {
            BucketSize[p-1]++;
            idx[k] = p-1;
        }
        else
        {
            for (i = 1; i < p-1; i++)
            {
                if (S[i-1] < data[k] and S[i] >= data[k])
                {
                    BucketSize[i]++;
                    idx[k] = i;
                    break;
                }
            }
        }
    }
    
    // Split A into p partitions B_0, B_1, B_2, .., B_p-1
    // Here Index[i] is used to get the current working index of the ith Bucket
    uint32_t Index[p];
    for (i = 0; i < p; i++)
    {
        Index[i] = 0;
    }

    Partition *B[p];
    for (i = 0; i < p; i++)
    {
        B[i] = new Partition(BucketSize[i]);
    }

    k = 0;
    for (k = 0; k < n; k++)
    {
        i = idx[k];
        B[i]->array[Index[i]] = data[k];
        Index[i]++;
    }
    
    for (i = 0; i < p; i++)
    {
        #pragma omp task
        {
            ParallelSort(B[i]->array, BucketSize[i], p);
        }
    }
    #pragma omp taskwait

    uint32_t *OffsetSize =  new uint32_t[p];
    for (i = 0; i < p; i++)
    {
        OffsetSize[i] = 0;
    }
    uint32_t offset = 0;

    for (i = 1 ; i < p; i++)
    {
        OffsetSize[i] = offset + BucketSize[i-1];
        offset = OffsetSize[i];
    }   

    k = 0;
    uint32_t r = 0;
    for (i = 0; i < p; i++)
    {
        #pragma omp task firstprivate(k)
        
        for (r = 0; r < BucketSize[i]; r++)
        {
            data[k+OffsetSize[i]] = B[i]->array[r];
            k++;
        }
    }
    #pragma omp taskwait
}
