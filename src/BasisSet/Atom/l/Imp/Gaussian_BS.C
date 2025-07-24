// File: Atom/l/Gaussian_BS.H
module;
#include <cmath>
#include <vector>
#include <iostream>
module qchem.BasisSet.Atom.Internal.l.GaussianBS;
import qchem.BasisSet.Atom.radial.Gaussian.ExponentScaler; 

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
