#include "enum_integer_points.hpp"



namespace SU2_Compiler{

/**
 * @brief Performs Gram–Schmidt orthogonalization on the columns of a basis matrix.
 *
 * Given an input basis matrix \p B whose columns are the original (possibly non-orthogonal)
 * basis vectors, this function returns a matrix of the same dimensions whose columns
 * form an orthonormal basis computed via the Gram–Schmidt process.
 *
 * @param[in] B  
 *     An m×n matrix whose n columns are the input basis vectors in ℝᵐ.
 * @return 
 *     An m×n matrix whose columns are the orthogonalized vectors corresponding to \p B.
 *     Note: columns are not scaled to unit length.
 *
 * @see  
 *     https://en.wikipedia.org/wiki/Gram–Schmidt_process
 */
MatrixXf GSO(const MatrixXf& B)
{
    int M = B.rows();
    int N = B.cols();
    MatrixXf B_orth(M, N);
    auto mu = [&](int i, int j) {
        return B.col(i).dot(B_orth.col(j)) / B_orth.col(j).squaredNorm();
    }; 

    for(int i = 0; i < N; i++){
        B_orth.col(i) = B.col(i);
        for(int j = 0; j < i; j++){
            B_orth.col(i) -= mu(i, j) * B_orth.col(j);
        }
    }
    return B_orth;
}



MatrixXi LLL(MatrixXf B, FTYPE delta)
{
    int M = B.rows();
    int N = B.cols();
    MatrixXf B_orth(M, N);
    MatrixXf mu(N, N);
    VectorXf norm2_B_orht(N);

    B_orth = GSO(B);
    for(int i = 0; i < N; i++) norm2_B_orht[i] = B_orth.squaredNorm();
    for(int i = 0; i < N; i++){
        for(int j = 0; j < i; j++){
            mu(i, j) = B.col(i).dot(B_orth.col(j)) / B_orth.col(j).squaredNorm();
        }
    }

    MatrixXi U = MatrixXf::Identity(N, N);

    int k = 1;
    while(k < N){
        for(int j = k-1; j >= 0; j--){
            if(abs(mu(k,j)) > 0.5){
                ITYPE q = (ITYPE)round(mu(k,j));
                B.col(k) -= q * B.col(j);
                U.col(k) -= q * U.col(j);
                for(int l = 0; l <= j; l++) mu(k,l) -= q * mu(j,l);
            }
        }
    
        if(norm2_B_orth(k) >= (delta - mu(k,k-1)*mu(k,k-1)) * norm2_B_orth(k-1)){
            k++;
        }else{
            B.col(k-1).swap(B.col(k));
            U.col(k-1).swap(U.col(k));

            FTYPE mu_prime = mu(k,k-1);
            FTYPE new_norm2_B = norm2_B_orth(k) + mu_prime*mu_prime * norm2_B_orth(k-1);
            mu(k,k-1) = mu_prime * norm2_B_orht(k-1) / new_norm2_B;
            norm2_B_orht(k) = norm2_B_orth(k) * norm2_B_orth(k-1) / new_norm2_B;
            norm2_B_orth(k-1) = new_norm2_B;
        }
    }

    
}

}