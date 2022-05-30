#include "classify.h"
#include <omp.h>
#include <stdio.h>
#include <string.h>

Data classify(Data &D, const Ranges &R, unsigned int numt)
{ 

   assert(numt < MAXTHREADS);
   Counter counts[R.num()]; 

   int index[D.ndata];
   #pragma omp parallel num_threads(numt) 
   {
      int tid = omp_get_thread_num(); // I am thread number tid
      #pragma omp for schedule(static, 50)
      for(int i=0; i<D.ndata; i+=numt) { // Threads together share-loop through all of Data
         int temp[numt];
         for (int j = 0; j < numt && i+j < D.ndata; j++)
         {
            int v = D.data[i+j].value = R.range(D.data[i+j].key);// For each data, find the interval of data's key,
                       // and store the interval id in value. D is changed.
            counts[v].increase(tid); // Found one key in interval v
            temp[j] = v;
         }
         memcpy(index + i, temp, numt * sizeof(int));
      }
   }
   
   unsigned int *rangecount = new unsigned int[R.num()];
   #pragma omp parallel num_threads(numt)
   {
      #pragma omp for
      for(int r=0; r<R.num(); r++) { 
         rangecount[r] = 0;
         for(int t=0; t<numt; t++) 
            rangecount[r] += counts[r].get(t);
         //std::cout << rangecount[r] << " elements in Range " << r << "\n"; // Debugging statement
      }
   }


   // Compute prefx sum on rangecount.
   for(int i=1; i<R.num(); i++) {
      rangecount[i] += rangecount[i-1];
   }

   Data D2 = Data(D.ndata); // Make a copy
   int rcount[R.num()] = {0};
   #pragma omp parallel num_threads(numt)
   {
      #pragma omp for
      for(int d=0; d<D.ndata; d++) // For each interval, thread loops through all of data and  
      {
         int idx = index[d];
         D2.data[rangecount[idx-1]+rcount[idx]] = D.data[d]; // Copy it to the appropriate place in D2.
         rcount[idx]++;
      }
   }

   return D2;
}
