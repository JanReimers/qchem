// File: Atom/l/BSpline_BS.H BSpline Basis Set for atoms.

#include "Imp/BasisSet/Atom/l/BSpline_BS.H"
#include "Imp/BasisSet/Atom/l/BSpline_IBS.H"
//#include "Imp/BasisSet/Atom/radial/BSpline/ExponentScaler.H"
//#include "Imp/BasisSet/Atom/radial/BSpline/Rk.H"

namespace Atoml
{
namespace BSpline
{


template <size_t K> BasisSet<K>::BasisSet(size_t N, double rmin, double rmax, size_t LMax)
{
    for (size_t L=0;L<=LMax;L++)
        Insert(new Orbital_IBS<K>(this,N,rmin,rmax,L));
        
}

template class BasisSet<6>;

}} //namespace
