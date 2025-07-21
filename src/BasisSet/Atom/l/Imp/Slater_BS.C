// File: Atom/l/Slater_BS.C Slater Basis Set for atoms.
module;
#include <vector>
#include <memory>
module qchem.BasisSet.Atom.l.SlaterBS;
import qchem.BasisSet.Atom.radial.Slater.ExponentScaler;
import qchem.BasisSet.Atom.radial.Slater.Rk;

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
