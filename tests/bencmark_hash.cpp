#include <iostream>
#include <unordered_map>
#include <vector>
#include <chrono>
#include <omp.h>
#include <execution>
#include <tbb/concurrent_unordered_map.h>
#include "tbb/concurrent_hash_map.h"

#include "src/HashTable.hpp"

using namespace std;
using namespace tbb;
using namespace SU2_Compiler;

using P = pair<long long, long long>;

template <>
struct std::hash<P>
{
    size_t operator()(const P& key) const
    {
        return key.first ^ (key.second << 1); 
    }
};

int main(){
    const int N = 1e8;
    vector<long long> counter(N);
    for(long long i = 0; i < N; i++) counter[i] = i;
    
    unordered_multimap<P, int> map_std;
    // concurrent_unordered_multimap<P, int> map_tbb;
    concurrent_hash_map<P, int> map_tbb(N);
    HashTable<P, int> my_hash(N);

    chrono::system_clock::time_point start, end;
    double time;

    const long long shift = 14343535;
    
    cout << "スレッド数 " << omp_get_max_threads() << endl;

    // start = chrono::system_clock::now();
    // for(auto &i : counter){
    //         // cout << h(i) << endl;
    //         map_std.insert({{i, i*i}, i});
    // }
    // end = chrono::system_clock::now();
    // time = static_cast<double>(chrono::duration_cast<chrono::microseconds>(end - start).count() / 1000.0);
    // cout << "std::mapの挿入時間 " << time << "[ms]" << endl;


    start = chrono::system_clock::now();
    tbb::parallel_for(0, N, [&](int i) {
        map_tbb.insert({{i, i*i}, i});
    });
// #pragma omp parallel for
//     for(auto &i : counter){
//             map_tbb.insert({{i, i*i}, i});
//     }
    end = chrono::system_clock::now();
    time = static_cast<double>(chrono::duration_cast<chrono::microseconds>(end - start).count() / 1000.0);
    cout << "tbb::mapの挿入時間 " << time << "[ms]" << endl;


    start = chrono::system_clock::now();
#pragma omp parallel for
    for(auto &i : counter){
            my_hash.insert({i, i*i}, i);
    }
    end = chrono::system_clock::now();
    time = static_cast<double>(chrono::duration_cast<chrono::microseconds>(end - start).count() / 1000.0);
    cout << "自作のmapの挿入時間 " << time << "[ms]" << endl;

}