// File: Atom/l/Gaussian_BS.H
module;
#include <cmath>
#include <vector>
#include <iostream>
#include "radial/Gaussian/ExponentScaler.H"
module qchem.BasisSet.Atom.l.GaussianBS;

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
