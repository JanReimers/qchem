// File: Atom/l/Slater_BS.C Slater Basis Set for atoms.

#include <vector>
#include "l/Slater_BS.H"
#include "l/Slater_IBS.H"
#include "radial/Slater/ExponentScaler.H"
#include "radial/Slater/Rk.H"

namespace Atoml
{
namespace Slater
{


BasisSet::BasisSet(size_t N, double emin, double emax, size_t LMax)
{
    ::Slater::ExponentScaler ss(N,emin,emax,LMax);
    for (size_t L=0;L<=LMax;L++)
        Insert(new Orbital_IBS(this,ss.Get_es(L),L));
        
}



}} //namespace
