#ifndef  CUDA_SUM
#define CUDA_SUM

#include<stdio.h>
#include<stdlib.h>
#include<time.h>

namespace plb {
namespace npfem {

template <class T>
__forceinline__ __device__ void sum(T *buffer, int n_sum, int n, int id) {


	for (int i = n_sum; i < n - id; i += n_sum) {
			buffer[id] += buffer[id + i];
	}
	int r = n_sum / 2;
	__syncthreads();

	while (r) {
		__syncthreads();
		if (id < r) {
			buffer[id] += buffer[id + r];
		}
		r /= 2;
	}
}

template <class T,  class T2>
__forceinline__ __device__ void double_sum(T *buffer, T2 *buffer2, int n_sum, int n, int id) {

	for (int i = n_sum; i < n - id; i += n_sum) {
		buffer[id] += buffer[id + i];
		buffer2[id] += buffer2[id + i];
	}
	int r = n_sum / 2;
	__syncthreads();
	while (r) {
		__syncthreads();
		if (id < r ) {
			buffer[id  ] += buffer[id + r ];

		}else if (id < 2*r) {
			buffer2[id - r] += buffer2[id];
		}
		r /= 2;
	}
	
}

template <class T>
__forceinline__ __device__ void triple_sum(T *buffer, T *buffer2,T *buffer3, int n_sum, int n, int id) {

	for (int i = n_sum; i < n - id; i += n_sum) {
			buffer[id]  += buffer[id + i];
			buffer2[id] += buffer2[id + i];
			buffer3[id] += buffer3[id + i];
	}
	int r = n_sum / 2;
	__syncthreads();
	if (id < r) {
		buffer[id] += buffer[id + r];
	}
	else if (id < 2 * r) {
		buffer2[id - r] += buffer2[id];
	}
	r /= 2;
	while (r) {
		__syncthreads();
		if (id < r) {
			buffer[id] += buffer[id + r];
		}else if (id < 2*r) {
			buffer2[id - r] += buffer2[id];
		}else if (id < 4*r) {
			buffer3[id - 2*r] += buffer3[id];
		}
		r /= 2;
	}
	__syncthreads();
	if (id == 0){
		buffer3[0] += buffer3[1];
	}
}

__global__ void test_sum(){

	__shared__ double buffer[516];
	buffer[threadIdx.x] = 1;
	buffer[threadIdx.x + 258] = 2;
	__syncthreads();
	double_sum(buffer, buffer + 258, 256, 258, threadIdx.x);
	__syncthreads();
	if(threadIdx.x == 0) printf("test sum %f %f",buffer[0], buffer[258]);
}

__global__ void test_sum2() {
	__shared__ double buffer1[258];
	__shared__ double buffer2[258];
	__shared__ double buffer3[258];
	buffer1[threadIdx.x] = threadIdx.x;
	buffer2[threadIdx.x] = threadIdx.x;
	buffer3[threadIdx.x] = threadIdx.x;
	__syncthreads();
	triple_sum(buffer1, buffer2, buffer3, 256, 258, threadIdx.x);
	__syncthreads();
	if (threadIdx.x == 0) printf("test sum %f %f %f", buffer1[0], buffer2[0], buffer3[0]);
}

}
}

#endif

