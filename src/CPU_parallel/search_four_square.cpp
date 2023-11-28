#pragma once

#include <bits/stdc++.h>
#include <omp.h>

template<class RandomAccessIterator>
using Pair_Iterator = std::pair<RandomAccessIterator, RandomAccessIterator>;

template<class RandomAccessIterator>
using Tuple_Iterator = std::tuple<RandomAccessIterator, RandomAccessIterator, RandomAccessIterator, RandomAccessIterator>;

/*
x^2 + y^2 = keyとなる要素を探索する. 
*/
template<class RandomAccessIterator, typename T>
Pair_Iterator<RandomAccessIterator> search_two_square(
                            RandomAccessIterator first_X, RandomAccessIterator last_X,
                            RandomAccessIterator first_Y, RandomAccessIterator last_Y, 
                            const T& key)
{
    for(auto it1 = first_X; it1 != last_X; it1++){
        for(auto it2 = first_Y; it2 != last_Y; it2++){
            if((*it1)*(*it1)+(*it2)*(*it2) == key){
                return {it1, it2};
            }
        }
    }

    return {last_X, last_Y};
}


/*
ソート済み配列A,Bからa_i + b_j = keyとなるi,jを探索する.
*/
template<class RandomAccessIterator, typename T>
std::vector<Pair_Iterator<RandomAccessIterator>> two_points_technique(
                            RandomAccessIterator first1, RandomAccessIterator last1,
                            RandomAccessIterator first2, RandomAccessIterator last2,
                            const T& key)
{
    std::vector<Pair_Iterator<RandomAccessIterator>> ret;
    static const T zero = 0;
    auto left = first1;
    auto right = last2 - 1;
    for(left = first1; left != last1; left++){
        while(right != first2-1){
            T diff = key - *left - *right;
            if(diff == zero){
                ret.push_back({left, right});
            }
            if(diff < zero) right--;
            else break;
        }
    }

    return ret;
}



/*
x^2 + y^2 + z^2 + w^2 = keyとなる要素を探索する. 
そのような要素があるときはそのiteratorを返す. 
存在しなときは最後のiteratorを返す. 
*/
template<class RandomAccessIterator, typename T>
std::vector<Tuple_Iterator<RandomAccessIterator>> search_four_square(
                            RandomAccessIterator first_X, RandomAccessIterator last_X,
                            RandomAccessIterator first_Y, RandomAccessIterator last_Y,
                            RandomAccessIterator first_Z, RandomAccessIterator last_Z,
                            RandomAccessIterator first_W, RandomAccessIterator last_W,
                            const T& key)
{
    unsigned long N1 = std::distance(first_X, last_X);
    unsigned long N2 = std::distance(first_Y, last_Y);
    unsigned long N3 = std::distance(first_Z, last_Z);
    unsigned long N4 = std::distance(first_W, last_W);


    int i = 0;
    std::vector<T> X_squared(N1);
    for(auto it = first_X; it != last_X; it++){
        X_squared[i] = (*it)*(*it);
        i++;
    }

    i = 0;
    std::vector<T> Y_squared(N2);
    for(auto it = first_Y; it != last_Y; it++){
        Y_squared[i] = (*it)*(*it);
        i++;
    }

    i = 0;
    std::vector<T> Z_squared(N3);
    for(auto it = first_Z; it != last_Z; it++){
        Z_squared[i] = (*it)*(*it);
        i++;
    }

    i = 0;
    std::vector<T> W_squared(N4);
    for(auto it = first_W; it != last_W; it++){
        W_squared[i] = (*it)*(*it);
        i++;
    }

    // 事前にソートしておくと交換回数が減って後のソートが少しだけ高速になる
    std::sort(X_squared.begin(), X_squared.end());
    std::sort(Y_squared.begin(), Y_squared.end());
    std::sort(Z_squared.begin(), Z_squared.end());
    std::sort(W_squared.begin(), W_squared.end());

    const int split_count = omp_get_max_threads();   // 時間計測するとこれが一番速そう?
    const int interval_X = std::max(int(std::ceil(double(N1) / split_count)), 1);
    const int interval_Z = std::max(int(std::ceil(double(N3) / split_count)), 1);

    std::vector<T*> XY_list(split_count, NULL);
    std::vector<T*> ZW_list(split_count, NULL);

#pragma omp parallel for
    for(int nX = 0; nX < split_count; nX++){
        const int XY_size = (std::min(N1, interval_X*(nX+1)) - interval_X*nX) * N2;
        if(XY_size <= 0) continue;
        T *XY = new T[XY_size];

        for(int i = 0; i < std::min(N1 - nX*interval_X, interval_X); i++){
            for(int j = 0; j < N2; j++) XY[i*N2 + j] = X2[nX*interval_X + i] + Y2[j];
        } 

        // sort_for_Direct_Product_noPrallel(XY, XY + XY_size, std::min(interval_X, N1-interval_X*nX), N2);
        std::sort(XY, XY + XY_size);

        XY_list[nX] = XY;
    }

#pragma omp parallel for
    for(int nZ = 0; nZ < split_count; nZ++){
        const int ZW_size = (std::min(N3, interval_Z*(nZ+1)) - interval_Z*nZ) * N4;
        if(ZW_size <= 0) continue;
        T *ZW = new T[ZW_size];

        for(int i = 0; i < min(n3 - nZ*interval_Z, interval_Z); i++){
            for(int j = 0; j < n4; j++) ZW[i*n4 + j] = Z2[nZ*interval_Z + i] + W2[j];
        } 
        
        // sort_for_Direct_Product_noPrallel(ZW, ZW + ZW_size, std::min(interval_Z, N3-interval_Z*nZ), N4);
        std::sort(ZW, ZW + ZW_size);

        ZW_list[nZ] = ZW;
    }


    
    std::vector<std::vector<Pair_Iterator<T*>>> 
                    results(split_count, std::vector<Pair_Iterator<T*>>(split_count));

#pragma omp parallel for 
    for(int loop = 0; loop < split_count*split_count; loop++){
        int nX = loop / split_count;
        int nZ = loop % split_count;

        T *XY = XY_list[nX];
        T *ZW = ZW_list[nZ];

        std::vector<Pair_Iterator<T*>> solutions = 
                                two_points_technique(XY, XY + N1*N2, ZW, ZW + N3*N4, key);
        
        results[nX][nZ] = solutions;
    }

    std::vector<Pair_Iterator<T*>> solutions_total;
    for(int nX = 0; nX < split_count; nX++){
        for(int nZ = 0; nZ < split_count; nZ++){
            for(auto ele : results[nX][nZ]) solutions_total.push_back(ele); 
        }
    }

    std::vector<Tuple_Iterator<RandomAccessIterator>> ret;
    for(auto P : solutions_total){
        T* left = P.first;
        T* right = P.second;
        
        auto [x_ans, y_ans] = search_two_square(first_X, last_X, first_Y, last_Y, *left);
        auto [z_ans, w_ans] = search_two_square(first_Z, last_Z, first_W, last_W, *right);

        ret.push_back({x_ans, y_ans, z_ans, w_ans});
    }  

    for(int nX = 0; nX < split_count; nX++) delete[] XY_list[nX];
    for(int nZ = 0; nZ < split_count; nZ++) delete[] XY_list[nZ];

    return ret;
}