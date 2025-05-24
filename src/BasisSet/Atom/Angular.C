// File: Angular_Sym.H Common interface for various atomic (spherical) symmetries.

#include "Imp/BasisSet/Atom/Angular.H"
#include <Spin.H>   

int ElCounts_l::GetN(Spin s) const
{
    return s==Spin::Up ? (n+nu)/2 : (n-nu)/2;  
}