#include <complex>
#include <vector>
#include <Eigen/Core>


/*
四元数クラス(a + bi + cj + dk)

[[ a+ib, c+di],
 [-c+di, a-bi]]
四元数a+bi+cj+dkから上のような2次元複素正方行列への写像を考えるとこれは単射環準同型となる(ハミルトン積で)。
また、単位四元数はSU(2)に対応する。
*/
template <typename T>
class quaternion{
private:
    T a;
    T b;
    T c; 
    T d;
    
public:
    quaternion(T a, T b, T c, T d) : a(a), b(b), c(c), d(d) {}
    
    // ハミルトン積
    quaternion operator*(const quaternion<T>& o){
        return quaternion(a*o.a - b*o.b - c*o.c - d*o.d, 
                          a*o.b + b*o.a + c*o.d - d*o.c,
                          a*o.c - b*o.d + c*o.a + d*o.b,
                          a*o.d + b*o.c - c*o.b + d*o.a);
    }
    
    // 共役作用(ユニタリの随伴に相当)
    quaternion<T> conj() {return {a, -b, -c, -d}; }

    T norm() {return std::sqrt(a*a + b*b + c*c + d*d);}

    std::tuple<T, T, T, T> get_abcd() {return {a,b,c,d};}
    
    Eigen::Matrix<T, 4, 1> get_abcd_by_eigen(){
        Eigen::Matrix<T, 4, 1> ret;
        ret << a, b, c, d;
        return ret;
    }

    Eigen::Matrix<std::complex<T>, 2, 2> get_Matrix(){
        Eigen::Matrix<std::complex<T>, 2, 2> ret;
        std::complex<T> u(a, b);
        std::complex<T> v(-c, d);
        ret << u, -std::conj(v), 
               v,  std::conj(u);
        return ret;
    }

    // ユニタリのトレース
    T trace() {return {2*a};}

    // 四元数がユニタリ(ノルムが1)かをチェック
    bool is_unitary(T tol=1e-14){
        if(abs(this->norm() - (T)1.0) < tol) return true;
        else return false;    
    }

    friend std::ostream& operator<<(std::ostream& os, const quaternion<T>& u){
        os << "(" << u.a << ", " << u.b << ", " << u.c << ", " << u.d << ")";
        return os;       
    }
};

template<typename T>
quaternion<T> operator*(const T& x, const quaternion<T>& y){
    return quaternion(x*y.a, x*y.b, x*y.c, x*y.d);
}

template<typename T>
quaternion<T> operator*(const quaternion<T>& x, const T& y){
    return quaternion(x.a*y, x.b*y, x.c*y, x.d*y);
}


template<typename T>
T distance(quaternion<T> u, quaternion<T> v){
    quaternion<T> uv_dag = u * v.conj();
    T tr = uv_dag.trace();
    return std::sqrt(4 - tr*tr);
}

template<typename T>
quaternion<T> convert_quaternion(Eigen::Matrix<T, 4, 1> v){
    quaternion<T> ret(v(0), v(1), v(2), v(3));
    return ret;
} 
