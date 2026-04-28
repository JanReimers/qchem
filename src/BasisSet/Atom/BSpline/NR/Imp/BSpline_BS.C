// File: Atom/l/BSpline_BS.C BSpline Basis Set for atoms.
module;
module qchem.BasisSet.Atom.BSpline.NR.BS;
import qchem.Symmetry.AtomEC;

namespace AtomBS
{
namespace BSpline
{
template <size_t K> BasisSet<K>::BasisSet(size_t N, double rmin, double rmax, const ElectronConfiguration& ec)
: AtomIE_BS_HF(this)
{
    const Atom_EC& aec=dynamic_cast<const Atom_EC&>(ec);
    size_t LMax=aec.GetLMax();
    for (auto ir:aec.GetIrreps())
        Insert(new Orbital_IBS<K>(this,N,rmin,rmax,ir));  
     
    BuildCache(LMax);
}

template <size_t K> void BasisSet<K>::Insert(Orbital_IBS<K>* oibs)
{
    ::BS_Common::Insert(oibs);
    AtomIE_BS_HF<double>::Append(oibs,oibs); //implicit casts to two different intefraces.
}

template <size_t K> BasisSet_r<K>::BasisSet_r(size_t N, double rmin, double rmax, const ElectronConfiguration& ec)
: AtomIE_BS_HF(this)
{
    const Atom_EC& aec=dynamic_cast<const Atom_EC&>(ec);
    size_t LMax=aec.GetLMax();
    for (auto ir:aec.GetIrreps())
         Insert(new Orbital_IBS_r<K>(this,N,rmin,rmax,ir));  
    BuildCache(LMax);
}

template <size_t K> void BasisSet_r<K>::Insert(Orbital_IBS_r<K>* oibs)
{
    ::BS_Common::Insert(oibs);
    AtomIE_BS_HF<double>::Append(oibs,oibs); //implicit casts to two different intefraces.
}

#define INSTANCEk(K) template class BasisSet<K>;
#include "../../Instance.hpp"
#define INSTANCEk(K) template class BasisSet_r<K>;
#include "../../Instance.hpp"


}} //namespace

