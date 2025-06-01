// File:  ElectronCounts.C  Simple structure that store electron configuration counts for each l state.

#include "Imp/BasisSet/Atom/ElectronCounts.H"
#include <Spin.H>
#include <cassert>


int ElCounts_l::GetN(Spin s) const
{
    return s==Spin::Up ? (N+Nu)/2 : (N-Nu)/2;  
}

void ElCounts::DebugCheck() const
{
    for (size_t l=0;l<=LMax;l++)
    {
        assert(Nv[l]>=Nu[l]);
        assert((Nv[l]-Nu[l])%2==0); 
        assert(N[l]==Nv[l]+Nf[l]);
        #ifndef NDEBUG
        size_t g=(2*l+1); //degenracy
        #endif
        assert(g>=Nu[l]);
        assert(2*g>=Nv[l]);
    }
}
    