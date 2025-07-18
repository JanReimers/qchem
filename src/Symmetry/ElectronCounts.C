// File:  ElectronCounts.C  Simple structure that store electron configuration counts for each l state.

#include <cassert>
#include <cstddef>
#include "Symmetry/ElectronCounts.H"
import qchem.Symmetry.Spin;


int ElCounts_l::GetN(Spin s) const
{
    assert(s!=Spin::None);
    return s==Spin::Up ? (N+Nu)/2 : (N-Nu)/2;  
}

int ElCounts::GetNv(int l, Spin s) const
{
    assert(s!=Spin::None);
    return s==Spin::Up ? (Nv[l]+Nu[l])/2 : (Nv[l]-Nu[l])/2;  
}

void ElCounts::DebugCheck() const
{
    for (size_t l=0;l<=LMax;l++)
    {
        assert(Nv[l]>=Nu[l]);
        assert((Nv[l]-Nu[l])%2==0); 
        assert(N[l]==Nv[l]+Nf[l]);
        #ifndef NDEBUG
        int g=(2*l+1); //degenracy
        #endif
        assert(g>=Nu[l]);
        assert(2*g>=Nv[l]);
        assert(Nf[l]%(2*g)==0);
    }
}
    