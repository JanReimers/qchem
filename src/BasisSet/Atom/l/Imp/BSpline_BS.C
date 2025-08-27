// File: Atom/l/BSpline_BS.C BSpline Basis Set for atoms.
module;
module qchem.BasisSet.Atom.Internal.l.BSplineBS;
import qchem.Symmetry.AtomEC;

namespace Atoml
{
namespace BSpline
{
template <size_t K> BasisSet<K>::BasisSet(size_t N, double rmin, double rmax, size_t LMax)
: Atom::BS_Common(this)
{
    for (size_t L=0;L<=LMax;L++)
        Insert(new Orbital_IBS<K>(this,N,rmin,rmax,L));
    BuildCache(LMax);
}


template <size_t K> BasisSet<K>::BasisSet(size_t N, double rmin, double rmax, const ElectronConfiguration& ec)
: Atom::BS_Common(this)
{
    const Atom_EC& aec=dynamic_cast<const Atom_EC&>(ec);
    size_t LMax=aec.GetLMax();
    for (size_t L=0;L<=LMax;L++)
    {
        auto mls=aec.GetBreadown(L);
        if (mls.ml_paired.size()>0)   
            Insert(new Orbital_IBS<K>(this,N,rmin,rmax,L,mls.ml_paired));            
        if (mls.ml_unpaired.size()>0)   
            Insert(new Orbital_IBS<K>(this,N,rmin,rmax,L,mls.ml_unpaired));            
        if (mls.ml_unoccupied.size()>0)   
            Insert(new Orbital_IBS<K>(this,N,rmin,rmax,L,mls.ml_unoccupied));            

    
    }           
    BuildCache(LMax);
}


#define INSTANCEk(k) template class BasisSet<k>;
#include "../../radial/BSpline/Instance.hpp"


}} //namespace

