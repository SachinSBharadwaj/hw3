#include <algorithm>
#include <stdio.h>
#include <math.h>
#include <omp.h>


// Scan A array and write result into prefix_sum array;
// use long data type to avoid overflow
void scan_seq(long* prefix_sum, const long* A, long n) {
  if (n == 0) return;
  prefix_sum[0] = A[0];
  for (long i = 1; i < n; i++) {
    prefix_sum[i] = prefix_sum[i-1] + A[i];
  }
	
}

void scan_omp(long* prefix_sum, const long* A, long n) {
  // TODO: implement multi-threaded OpenMP scan
	if (n == 0) return;
	prefix_sum[0] = A[0];
	long NT = 8; // # of threads
	long jump = n/NT;
	long* last = (long*) malloc(NT * sizeof(long));
	long* lsum = (long*) malloc(NT * sizeof(long));
	last[0]=0;
	long i,k;

	#pragma omp parallel for num_threads(NT) 
	for (long k = 0; k < NT; k = k+1){
		prefix_sum[k*jump]=A[k*jump];
		#pragma omp parallel for num_threads(2) 
		for (long i = 1; i < jump; i++) {
			
    			prefix_sum[k*jump+i] = prefix_sum[k*jump+i-1] + A[k*jump+i];
			if(i==(jump-1) && k!=(NT-1)){last[k+1]=prefix_sum[k*jump+i];}
			
  		}
	}
	
	#pragma omp barrier
	
	scan_seq(lsum, last, NT);

	#pragma omp parallel for num_threads(NT)
	for (long k = 0; k < NT;k = k+1){
		for (long i =0; i<jump; i=i+1){
			prefix_sum[k*jump+i] = prefix_sum[k*jump+i] + lsum[k];		
		}	
	}
	
}

int main() {
  long N = 100000000;
  long* A = (long*) malloc(N * sizeof(long));
  long* B0 = (long*) malloc(N * sizeof(long));
  long* B1 = (long*) malloc(N * sizeof(long));
  for (long i = 0; i < N; i++) {A[i] = rand();}

  double tt = omp_get_wtime();
  scan_seq(B0, A, N);
  printf("sequential-scan = %fs\n", omp_get_wtime() - tt);

  tt = omp_get_wtime();
  scan_omp(B1, A, N);
  printf("parallel-scan   = %fs\n", omp_get_wtime() - tt);

  long err = 0;
  for (long i = 0; i < N; i++) err = std::max(err, std::abs(B0[i] - B1[i]));
  printf("error = %ld\n", err);

  free(A);
  free(B0);
  free(B1);
  return 0;
}
