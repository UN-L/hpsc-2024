#include <cstdio>
#include <cstdlib>
#include <vector>

__global__ void bucketSort(int range, int bucket[], int key[], int offset[])
{
  int i = threadIdx.x;
  for (int j = 1; j<range; j<<=1) {
    offset[i] = bucket[i];
    if(i>=j) bucket[i] += offset[i-j];
  }

  for (int j=0; bucket[i]>0; bucket[i]--)
    key[j++] = i;
}

int main() 
{
  int n = 50;
  int range = 5;
  int *key;
  cudaMallocManaged(&key, n*sizeof(int));

  for (int i=0; i<n; i++) {
    key[i] = rand() % range;
    printf("%d ",key[i]);
  }

  printf("\n");

  int *bucket;
  cudaMallocManaged(&bucket, range*sizeof(int));
  int *offset;
  cudaMallocManaged(&offset, range*sizeof(int)); 
  for (int i=0; i<range; i++) {
    bucket[i] = 0;
  }
  for (int i=0; i<n; i++) {
    bucket[key[i]]++;
  }
  bucketSort<<<1, n>>>(range, bucket, key, offset);
  cudaDeviceSynchronize();

  for (int i=0; i<n; i++) {
    printf("%d ",key[i]);
  }
  printf("\n");
  cudaFree(key);
  cudaFree(bucket);
}
