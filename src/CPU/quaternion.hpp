#pragma once

#include <complex>
#include <vector>
#include <Eigen/Core>
#include "Rings.cpp"


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

    inline T get_a() const {return a;}
    inline T get_b() const {return b;}
    inline T get_c() const {return c;}
    inline T get_d() const {return d;}
    
    // ハミルトン積
    inline quaternion operator*(const quaternion<T>& o){
        return quaternion(a*o.a - b*o.b - c*o.c - d*o.d, 
                          a*o.b + b*o.a + c*o.d - d*o.c,
                          a*o.c - b*o.d + c*o.a + d*o.b,
                          a*o.d + b*o.c - c*o.b + d*o.a);
    }

    inline quaternion operator+(const quaternion<T>& o){
        return quaternion(a+o.a, b+o.b, c+o.c, d+o.d);
    }

    inline quaternion operator-(const quaternion<T>& o){
        return quaternion(a-o.a, b-o.b, c-o.c, d-o.d);
    }


    inline quaternion operator-(){
        return quaternion(-a,-b,-c,-d);
    }

    inline quaternion operator/(const T& o){
        return quaternion(a / o, b / o, c / o, d / o);
    }
    
    // 共役作用(ユニタリの随伴に相当)
    quaternion<T> conj() {return {a, -b, -c, -d}; }

    T norm() {return sqrt(a*a + b*b + c*c + d*d);}

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

    void unitalize(){
        T n = this->norm();
        a /= n;
        b /= n;
        c /= n;
        d /= n;
    }

    friend std::ostream& operator<<(std::ostream& os, const quaternion<T>& u){
        os << "(" << u.a << ", " << u.b << ", " << u.c << ", " << u.d << ")";
        return os;       
    }
};

template<typename T>
inline quaternion<T> operator*(const T& x, const quaternion<T>& y){
    return quaternion(x*y.get_a(), x*y.get_b(), x*y.get_c(), x*y.get_d());
}

template<typename T>
inline quaternion<T> operator*(const quaternion<T>& x, const T& y){
    return quaternion(x.get_a()*y, x.get_b()*y, x.get_c()*y, x.get_d()*y);
}


// 1/2||u - v||◇を計算
template<typename T>
T distance(quaternion<T> u, quaternion<T> v){
    quaternion<T> uv_dag = u * v.conj();
    T tr = uv_dag.trace();
    return sqrt(1 - (tr*tr / 4.0));
}

// u-vの最大値ノルムを計算
template <typename T>
T distance_max(quaternion<T> u, quaternion<T> v){
    quaternion<T> diff = u -v;
    std::array<T, 4> arr = {abs(diff.get_a()), abs(diff.get_b()), abs(diff.get_c()), abs(diff.get_d())};
    return *std::max_element(arr.begin(), arr.end());
}

template<typename T>
quaternion<T> convert_quaternion(Eigen::Matrix<T, 4, 1> v){
    quaternion<T> ret(v(0), v(1), v(2), v(3));
    return ret;
} 


template<typename ITYPE, typename FTYPE>
quaternion<FTYPE> to_quaterion(ZOmega<ITYPE> x, ZOmega<ITYPE> y){
    const FTYPE sqrt2 = sqrt((FTYPE)2.0);
    FTYPE a = (FTYPE)x.d + (FTYPE)(x.c - x.a) * sqrt2 / 2.0;
    FTYPE b = (FTYPE)x.b + (FTYPE)(x.c + x.a) * sqrt2 / 2.0;
    FTYPE c = (FTYPE)y.d + (FTYPE)(y.c - y.a) * sqrt2 / 2.0;
    FTYPE d = (FTYPE)y.b + (FTYPE)(y.c + y.a) * sqrt2 / 2.0;
    return quaternion<FTYPE>(a,b,c,d);
}


template<typename T>
quaternion<T> random_unitary(size_t seed=-1){
    std::random_device rd;
    std::default_random_engine eng(rd());
    if(seed != -1) eng.seed(seed);
    std::uniform_real_distribution<double> distr(-1.0, 1.0);

    T a = distr(eng);
    T b = distr(eng);
    T c = distr(eng);
    T d = distr(eng);
    quaternion<T> runitary(a,b,c,d);
    runitary.unitalize();
    return runitary;
}