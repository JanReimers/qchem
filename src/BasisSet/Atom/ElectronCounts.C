// File:  ElectronCounts.C  Simple structure that store electron configuration counts for each l state.

#include "Imp/BasisSet/Atom/ElectronCounts.H"
#include <Spin.H>


int ElCounts_l::GetN(Spin s) const
{
    return s==Spin::Up ? (N+Nu)/2 : (N-Nu)/2;  
}