// File: IntPower.H  real**int functions since C doesn't provide them.
module;
#include <iostream>
#include <stdlib.h>

export module Common.IntPow;

export 
{
    inline double uintpow(double x, unsigned int n)
    {
        if (n==0) return 1.0;
        if (x==0) return 0.0;
        if (n==1) return x;
        if (n==2) return x*x;
        if (n==3) return x*x*x;
        if (n==4) return (x*x)*(x*x);
        double ret=uintpow(x,n/2);
        if (n%2)
        {
            return x*ret*ret;
        }
        else
        {
            return ret*ret;
        }

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