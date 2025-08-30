// File: Atom/l/BSpline_BS.C BSpline Basis Set for atoms.
module;
module qchem.BasisSet.Atom.BSpline.NR.BS;
import qchem.Symmetry.AtomEC;

namespace Atom
{
namespace BSpline
{
template <size_t K> BasisSet<K>::BasisSet(size_t N, double rmin, double rmax, size_t LMax)
: AtomIE_BS_2E(this)
{
    for (size_t L=0;L<=LMax;L++)
        Insert(new Orbital_IBS<K>(this,N,rmin,rmax,L));
    BuildCache(LMax);
}


template <size_t K> BasisSet<K>::BasisSet(size_t N, double rmin, double rmax, const ElectronConfiguration& ec)
: AtomIE_BS_2E(this)
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

template <size_t K> void BasisSet<K>::Insert(Orbital_IBS<K>* oibs)
{
    ::BS_Common::Insert(oibs);
    AtomIE_BS_2E<double>::Append(oibs,oibs); //implicit casts to two different intefraces.
}

#define INSTANCEk(k) template class BasisSet<k>;
#include "../../Instance.hpp"


}} //namespace

