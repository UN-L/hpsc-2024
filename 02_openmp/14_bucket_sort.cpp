#include <cstdio>
#include <cstdlib>
#include <vector>
#include <omp.h>

int main() {
  int n = 50000000;
  int range = 5;
  std::vector<int> key(n);
  for (int i=0; i<n; i++) {
    key[i] = rand() % range;
    printf("%d ",key[i]);
  }
  printf("\n");

  std::vector<int> bucket(range,0); 
  for (int i=0; i<n; i++)
    bucket[key[i]]++;
  std::vector<int> offset(range,0);
  for (int i=1; i<range; i++) 
    offset[i] = offset[i-1] + bucket[i-1];
  
  //Compile but not work
  #pragma omp parallel for
  for (int i=0; i<range; i++) {
    int j = offset[i];
    for (int k=bucket[i]; k>0; k--) {
      key[j++] = i;
    }
  }

  //I've tried separating for loops naively but I don't know how to separate offset in different threads 
  /*#pragma omp parallel 
  {
    int i =0;
    int j =offset[i];
    std::vector<int> tmp_vec(range,0);
    #pragma omp for 
    for (i=0; i<range; i++) {
      tmp_vec[i] = offset[i];
    }
    #pragma omp for 
    for (int k=bucket[i]; k>0; k--) {
      key[j++] = i;
    }
  }*/

  for (int i=0; i<n; i++) {
    printf("%d ",key[i]);
  }
  printf("\n");
}