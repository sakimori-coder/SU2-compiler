#include <iostream>
#include <algorithm>        // For STL sort
#include <thrust/sort.h>    // For thrust sort
#include <thrust/binary_search.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <time.h>

using namespace std;

struct my_type
{
    int a;
    int b;

    // __host__ __device__ my_type inline operator=(my_type other){
    //     a = other.a;
    //     b = other.b;
    //     return *this;
    // }
};


__host__ __device__ bool cmp(pair<int, int> x, pair<int, int> y){
    return x.first < y.first;
}

int main()
{
    // Creating vector with random numbers
    const int N = 100000000;
    thrust::host_vector<my_type> A_host(N), B_host(N);
    for (int i = 0; i < N; i++) {
        my_type tmp = {rand(), rand()};
        A_host[i] = tmp;
        B_host[i] = tmp;
    }

    // Data transfer (CPU -> GPU) 
    thrust::device_vector<my_type> A(N), B(N);
    thrust::copy(A_host.begin(), A_host.end(), A.begin());
    thrust::copy(B_host.begin(), B_host.end(), B.begin());


    clock_t start = clock();
    
    // Sort (thrust)
    // thrust::sort(A.begin(), A.end(), cmp);   // Sort by thrust
    // thrust::sort(B.begin(), B.end());
    clock_t end = clock();
    double time = static_cast<double>(end - start) / CLOCKS_PER_SEC * 1000.0;
    printf("time %lf[ms]\n", time);
    // Data transfer (GPU -> CPU)
    // thrust::copy(d_vec.begin(), d_vec.end(), h_vec.begin());

    return 0;
}