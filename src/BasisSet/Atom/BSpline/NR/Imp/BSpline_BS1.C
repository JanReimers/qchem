// File: Atom/l/BSpline_BS.C BSpline Basis Set for atoms.
module;
module qchem.BasisSet.Atom.BSpline.NR.BS1;
import qchem.Symmetry.AtomEC;

namespace AtomBS
{
namespace BSpline
{
template <size_t K> BasisSet1<K>::BasisSet1(size_t N, double rmin, double rmax, const ElectronConfiguration& ec)

{
    const Atom_EC& aec=dynamic_cast<const Atom_EC&>(ec);
    size_t LMax=aec.GetLMax();
    for (auto ir:aec.GetIrreps())
        Insert(new Orbital_IBS1<K>(N,rmin,rmax,ir));  
     
    BuildCache(LMax);
}

template <size_t K> void BasisSet1<K>::Insert(Orbital_IBS1<K>* oibs)
{
    ::BS_Common1::Insert(oibs);
    // AtomIE_BS_HF<double>::Append(oibs,oibs); //implicit casts to two different intefraces.
}

template <size_t K> BasisSet_r1<K>::BasisSet_r1(size_t N, double rmin, double rmax, const ElectronConfiguration& ec)
// : AtomIE_BS_HF(this)
{
    const Atom_EC& aec=dynamic_cast<const Atom_EC&>(ec);
    size_t LMax=aec.GetLMax();
    for (auto ir:aec.GetIrreps())
         Insert(new Orbital_IBS_r1<K>(N,rmin,rmax,ir));  
    BuildCache(LMax);
}

template <size_t K> void BasisSet_r1<K>::Insert(Orbital_IBS_r1<K>* oibs)
{
    ::BS_Common1::Insert(oibs);
    // AtomIE_BS_HF<double>::Append(oibs,oibs); //implicit casts to two different intefraces.
}

#define INSTANCEk(K) template class BasisSet1<K>;
#include "../../Instance.hpp"
#define INSTANCEk(K) template class BasisSet_r1<K>;
#include "../../Instance.hpp"


}} //namespace

