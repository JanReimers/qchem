// File Slater_m/BasisSet.H

#include "Imp/BasisSet/Atom/ml/Slater_BS.H"
#include "Imp/BasisSet/Atom/ml/Slater_IBS.H"
#include "Imp/BasisSet/Atom/radial/Slater/ExponentScaler.H"

namespace Atom_ml
{
namespace Slater
{


BasisSet::BasisSet(size_t N, double emin, double emax, size_t LMax)
{
    ::Slater::ExponentScaler ss(N,emin,emax,LMax);
    for (size_t L=0;L<=LMax;L++)
        for (int m=-L;m<=(int)L;m++)
            Insert(new Orbital_IBS(this,ss.Get_es(L),L,m));            

}

}} //namespace
