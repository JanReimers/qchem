// File: Atom/l/Slater_BS.C Slater Basis Set for atoms.
module;
#include <vector>
#include <memory>
#include "radial/Slater/ExponentScaler.H"
#include "radial/Slater/Rk.H"
module qchem.BasisSet.Atom.l.SlaterBS;
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
