#include "type.hpp"
#include "quaternion.hpp"

#include <random>

namespace su2compiler
{
    Matrix2cf quaternion::get_Matrix() const
    {
        Matrix2cf ret;
        CTYPE u(a,b);
        CTYPE t(c,d);
        ret << u, -std::conj(t), 
               t,  std::conj(u);
        
        return ret;
    }

    bool quaternion::is_unitary(FTYPE tol) const
    {
        if(abs(this->norm() - 1) < tol) return true;
        else return false;
    }

    void quaternion::unitalize()
    {
        FTYPE norm = this->norm();
        a /= norm;
        b /= norm;
        c /= norm;
        d /= norm;
    }

    // 通常のハミルトン積とは定義が少し違う
    quaternion& quaternion::operator*=(const quaternion& r)
    {
        FTYPE ret_a = a*r.a - b*r.b - c*r.c - d*r.d;
        FTYPE ret_b = a*r.b + b*r.a - c*r.d + d*r.c;
        FTYPE ret_c = a*r.c + b*r.d + c*r.a - d*r.b;
        FTYPE ret_d = a*r.d - b*r.c + c*r.b + d*r.a;
    
        a = ret_a;
        b = ret_b;
        c = ret_c;
        d = ret_d;

        return *this;
    }

    quaternion quaternion::operator*(const quaternion& r) const
    {
        quaternion ret = *this;
        ret *= r;
        return ret;
    }

    std::ostream& operator<<(std::ostream& os, const quaternion& U)
    {
        os << "(" << U.a << ", " << U.b << ", " << U.c << ", " << U.d << ")";
        return os;
    }

    quaternion adjoint(const quaternion& U)
    {
        return {U.a, -U.b, -U.c, -U.d};
    }

    FTYPE distance(const quaternion& U, const quaternion& V)
    {
        quaternion UV_dag = U * adjoint(V);
        FTYPE tr = UV_dag.trace();
        return sqrt(1.0 - (tr*tr) / 4.0);
    }

    std::random_device rd;
    std::default_random_engine eng(rd());
    
    void set_random_unitary_seed(size_t seed)
    {
        eng.seed(seed);
    }

    quaternion random_unitary()
    {
        std::uniform_real_distribution<double> distr(-1.0, 1.0);
        
        FTYPE a = distr(eng);
        FTYPE b = distr(eng);
        FTYPE c = distr(eng);
        FTYPE d = distr(eng);
        quaternion Rand_Unitary(a,b,c,d);
        Rand_Unitary.unitalize();
        return Rand_Unitary;
    }
}
