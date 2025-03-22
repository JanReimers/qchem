// File Slater_m/BasisSet.H

#include "Imp/BasisSet/Slater_m/BasisSet.H"
#include "Imp/BasisSet/Slater_m/IrrepBasisSet.H"
#include "Imp/BasisSet/SlaterScaler.H"

namespace Slater_m
{


BasisSet::BasisSet(size_t N, double emin, double emax, size_t LMax)
{
    SlaterScaler ss(N,emin,emax,LMax);
    for (size_t L=0;L<=LMax;L++)
        for (int m=-L;m<=(int)L;m++)
            Insert(new Orbital_IBS(this,ss.Get_es(L),L,m));            

}

} //namespace
