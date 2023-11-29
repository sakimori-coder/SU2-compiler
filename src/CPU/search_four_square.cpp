#pragma once

#include <bits/stdc++.h>

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

    T *XY, *ZW;
    XY = new T[N1*N2];
    ZW = new T[N3*N4];

    for(int i = 0; i < N1; i++){
        for(int j = 0; j < N2; j++) XY[i*N2 + j] = X_squared[i] + Y_squared[j];
    }
    for(int i = 0; i < N3; i++){
        for(int j = 0; j < N4; j++) ZW[i*N4 + j] = Z_squared[i] + W_squared[j];
    } 

    std::sort(XY, XY + N1*N2);
    std::sort(ZW, ZW + N3*N4);


    std::vector<Pair_Iterator<T*>> solutions = 
                            two_points_technique(XY, XY + N1*N2, ZW, ZW + N3*N4, key);


    std::vector<Tuple_Iterator<RandomAccessIterator>> ret;
    for(auto P : solutions){
        T* left = P.first;
        T* right = P.second;
        
        auto [x_ans, y_ans] = search_two_square(first_X, last_X, first_Y, last_Y, *left);
        auto [z_ans, w_ans] = search_two_square(first_Z, last_Z, first_W, last_W, *right);

        ret.push_back({x_ans, y_ans, z_ans, w_ans});
    }

    delete[] XY;
    delete[] ZW;

    return ret;
}




template<class RandomAccessIterator, typename T>
std::vector<Tuple_Iterator<RandomAccessIterator>> search_four_square_using_vector(
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

    std::vector<T> XY(N1*N2), ZW(N3*N4);

    for(int i = 0; i < N1; i++){
        for(int j = 0; j < N2; j++) XY[i*N2 + j] = X_squared[i] + Y_squared[j];
    }
    for(int i = 0; i < N3; i++){
        for(int j = 0; j < N4; j++) ZW[i*N4 + j] = Z_squared[i] + W_squared[j];
    } 

    std::sort(XY.begin(), XY.end());
    std::sort(ZW.begin(), ZW.end());


    auto solutions = two_points_technique(XY.begin(), XY.end(), ZW.begin(), ZW.end(), key);


    std::vector<Tuple_Iterator<RandomAccessIterator>> ret;
    for(auto P : solutions){
        auto left = P.first;
        auto right = P.second;
        
        auto [x_ans, y_ans] = search_two_square(first_X, last_X, first_Y, last_Y, *left);
        auto [z_ans, w_ans] = search_two_square(first_Z, last_Z, first_W, last_W, *right);

        ret.push_back({x_ans, y_ans, z_ans, w_ans});
    }

    return ret;
}