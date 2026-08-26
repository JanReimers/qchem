// File: IntPower.H  real**int functions since C doesn't provide them.
module;
#include <iostream>
#include <stdlib.h>

export module qchem.IntPow;

export 
{
    //! The n>4 tail of \c uintpow: binary powering, RECURSIVE.  Split out so that \c uintpow itself is
    //! NON-recursive and therefore INLINABLE -- clang will not inline a self-recursive function, and
    //! uintpow sits three deep in the collocation kernel's innermost loop (perf 2026-08-26: 7.3% of the
    //! uncached GPW box walk, half of it in the @plt stub, i.e. a cross-DSO call per Cartesian factor).
    //! Associations are UNCHANGED from the original single function, so every value is bit-identical.
    inline double uintpow_tail(double x, unsigned int n)
    {
        if (n==1) return x;
        if (n==2) return x*x;
        if (n==3) return x*x*x;
        if (n==4) return (x*x)*(x*x);
        double ret=uintpow_tail(x,n/2);
        if (n%2)
        {
            return x*ret*ret;
        }
        else
        {
            return ret*ret;
        }

    }

    inline double uintpow(double x, unsigned int n)
    {
        if (n==0) return 1.0;
        if (x==0) return 0.0;
        if (n==1) return x;
        if (n==2) return x*x;
        if (n==3) return x*x*x;
        if (n==4) return (x*x)*(x*x);
        return uintpow_tail(x,n);
    }

    inline double  intpow(double x,          int n)
    {
        if (n==0) return 1.0;
        if (x==0.0)
        {
            if (n>=0)
                return 0.0;
            else
            {
                std::cerr << "intpow::Divide by zero" << std::endl;
                exit(-1);
            }
        }
        double ret;
        if (n>0) ret=x;
        else ret=1.0/x;
        int np=abs(n);

        return uintpow(ret,np);
    }
} //export block