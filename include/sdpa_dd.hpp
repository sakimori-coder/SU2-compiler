#include <iostream>
#include <fstream>
#include <string>
#include <algorithm>
#include <cctype>
#include <Eigen/Core>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/algorithm/string.hpp>

namespace mp = boost::multiprecision;

template <typename T>
using vector_4d = std::vector<std::vector<std::vector<std::vector<T>>>>;


template<typename T>
void input_sdpa_data(std::vector<T>& C, vector_4d<T>& F, std::ofstream& ofs){
    ofs << std::setprecision(50);
    
    int mDIM = C.size();
    int nBLOCK = F[0].size();
    std::vector<int> bLOCKsTRUCT(nBLOCK);
    for(int i = 0; i < nBLOCK; i++){
        if(F[0][i].size() == F[0][i][0].size()) bLOCKsTRUCT[i] = F[0][i].size();
        else bLOCKsTRUCT[i] = -F[0][i][0].size();
    }

    ofs << mDIM << "  = mDIM" << std::endl;
    ofs << nBLOCK << " = nBLOCK" << std::endl;
    for(int block_size : bLOCKsTRUCT) ofs << block_size << " ";
    ofs << " = bLOCKsTRUCT" << std::endl;

    ofs << "{";
    for(int i = 0; i < mDIM; i++) ofs << C[i] << ", ";
    ofs << "}" << std::endl;

    
    for(int i = 0; i < mDIM+1; i++){
        ofs << "{";
        for(int j = 0; j < nBLOCK; j++){
            ofs << "{";
            if(bLOCKsTRUCT[j] < 0){
                ofs << "{";
                for(int k = 0; k < -bLOCKsTRUCT[j]; k++){
                    ofs << F[i][j][0][k] << ", ";
                }
                ofs << "}" << std::endl;
            }else{
                for(int k = 0; k < bLOCKsTRUCT[j]; k++){
                    ofs << "{"; 
                    for(int l = 0; l < bLOCKsTRUCT[j]; l++){
                        ofs << F[i][j][k][l] << ", ";
                    }
                    ofs << "}," << std::endl;
                }
            }
            ofs << "}" << std::endl;
        }
        ofs << "}" << std::endl;
    } 

    ofs << std::setprecision(6);
}

// string to float
template <typename T>
T stof(std::string s){
    T ret(s);
    return ret;
}

template<>
double stof(std::string s){
    return std::stod(s);
}


template <typename T>
std::vector<T> output_sdpa_data(std::ifstream& ifs){
    std::string line;
    while(getline(ifs, line)){
        if(line == "xVec = ") break;
    }

    getline(ifs, line);
    line = line.substr(1, line.size()-2);
    std::vector<std::string> xVec_str;
    boost::algorithm::split(xVec_str, line, boost::is_any_of(","));
    std::vector<T> xVec;
    for(std::string x_str : xVec_str){
        boost::algorithm::trim(x_str);
        T x = stof<T>(x_str);
        xVec.push_back(x);
    }   
    
    return xVec;
}
    

template <typename T>
std::vector<T> sdpa_solver(std::vector<T>& C, vector_4d<T>& F){
    std::string filename1 = "input_sdpa.dat";
    std::string filename2 = "output_sdpa.dat";
    std::string filename3 = "log_sdpa";
    std::ofstream ofs(filename1);
    if(ofs.fail()){
        std::cout << filename1 + "が開けません" << std::endl;
        exit(1);
    }

    input_sdpa_data(C, F, ofs);
    // std::string command_sdpa = "sdpa_gmp " + filename1 + " " + filename2 + " > " + filename3;
    std::string command_sdpa = "sdpa_qd " + filename1 + " " + filename2 + " > " + filename3; 
    if(system(command_sdpa.c_str()) == -1){
        std::cout << "SDPコマンドが正常に動作しませんでした" << std::endl;
    }

    std::ifstream ifs(filename2);
    if(ifs.fail()){
        std::cout << filename2 + "が開けません" << std::endl;
        exit(1);
    }
    std::vector<T> xVec = output_sdpa_data<T>(ifs);

    std::string command_rm = "rm " + filename1 + " " + filename2 + " " + filename3;
    // if(system(command_rm.c_str()) == -1){
    //     std::cout << filename1 << "と" << filename2 << "と" << filename3 << "が削除されていません" << std::endl;
    // }
    
    return xVec;
}