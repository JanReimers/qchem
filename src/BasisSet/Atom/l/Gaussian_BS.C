// File: Atom/l/Gaussian_BS.H

#include <cmath>
#include <vector>
#include <iostream>
#include "l/Gaussian_BS.H"
#include "l/Gaussian_IBS.H"
#include "radial/Gaussian/ExponentScaler.H"

namespace Atoml
{
namespace Gaussian
{

BasisSet::BasisSet(size_t N, double emin, double emax, size_t LMax)
{
    ::Gaussian::ExponentScaler gs(N,emin,emax,LMax);
    for (size_t L=0;L<=LMax;L++)
        Insert(new Orbital_IBS(this,gs.Get_es(L),L)); 
}



} //namespace
} //namespace
