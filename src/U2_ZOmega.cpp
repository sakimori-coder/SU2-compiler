#include "U2_ZOmega.hpp"
#include "rings.hpp"
#include "quaternion.hpp"
#include <iostream>


namespace SU2_Compiler
{
    Mat2ZOmega U2_ZOmega::get_Matrix() const
    {
        Mat2ZOmega ret;
        ret << u, -conj(t) * omega_pow[l],
               t,  conj(u) * omega_pow[l];
        
        return ret;
    }


    inline U2_ZOmega& U2_ZOmega::operator*=(const U2_ZOmega& r)
    {
        ZOmega omega = {0,0,1,0};
        ZOmega ret_u = u * r.u - conj(t) * r.t * omega;
        ZOmega ret_t = t * r.u + conj(u) * r.t * omega;
        int ret_l = (l + r.l) % 8;
        int ret_k = k + r.k;

        u = ret_u;
        t = ret_t;
        l = ret_l;
        k = ret_k;

        reduction(*this);
        
        return *this;
    }
    inline U2_ZOmega U2_ZOmega::operator*(const U2_ZOmega& r) const
    {
        U2_ZOmega ret = *this;
        ret *= r;
        return ret;
    }
    inline U2_ZOmega U2_ZOmega::operator-() const
    {
        return {u, t, l, k};
    }

    U2_ZOmega adjoint(const U2_ZOmega& U)
    {
        U2_ZOmega U_dag;
        U_dag.u = conj(U.u);
        U_dag.t = -U.t * conj(omega_pow[U.l]);
        U_dag.l = (8 - U.l) % 8;
        U_dag.k = U.k;
        return U_dag;
    }

    void reduction(U2_ZOmega& U)
    {
        while(true){
            ZOmega u = U.u;
            ZOmega t = U.t;
            bool flag_u = (u.a - u.c) % 2 || (u.b - u.d) % 2;
            bool flag_t = (t.a - t.c) % 2 || (t.b - t.d) % 2;
            if(flag_u || flag_t) return;

            U.k--;
            
            U.u.a = (u.b - u.d) / 2;
            U.u.b = (u.c + u.a) / 2;
            U.u.c = (u.b + u.d) / 2;
            U.u.d = (u.c - u.a) / 2;

            U.t.a = (t.b - t.d) / 2;
            U.t.b = (t.c + t.a) / 2;
            U.t.c = (t.b + t.d) / 2;
            U.t.d = (t.c - t.a) / 2;
        }
    }

    quaternion to_quaternion(const U2_ZOmega& U)
    {
        CTYPE phase = 1.0;
        for(int i = 0; i < U.l; i++) phase *= sqrt_omega;

        CTYPE u = ZOmega_to_FTYPE(U.u) * std::conj(phase);
        CTYPE t = ZOmega_to_FTYPE(U.t) * std::conj(phase);

        u /= pow(sqrt2, (FTYPE)U.k);
        t /= pow(sqrt2, (FTYPE)U.k);        
        
        quaternion ret(u.real(), u.imag(), t.real(), t.imag());
        return ret;
    }

    ZOmega omega_pow[8] = {{0,0,0, 1}, {0,0, 1,0}, {0, 1,0,0},  {1,0,0,0},
                           {0,0,0,-1}, {0,0,-1,0}, {0,-1,0,0}, {-1,0,0,0}};
}
