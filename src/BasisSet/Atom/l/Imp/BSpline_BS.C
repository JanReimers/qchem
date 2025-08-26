// File: Atom/l/BSpline_BS.C BSpline Basis Set for atoms.
module;
#include <iostream>
#include <cassert>

module qchem.BasisSet.Atom.Internal.l.BSplineBS;
import BasisSet.Atom.BSpline_IBS;
import qchem.BasisSet.Atom.Internal.radial.BSpline.Rk;
import qchem.BasisSet.Atom.Internal.Angular;
import qchem.Symmetry.AtomEC;

namespace Atoml
{
namespace BSpline
{
template <size_t K> BasisSet<K>::BasisSet(size_t N, double rmin, double rmax, size_t LMax)
: ::BSpline::BS_Common<K>(this,new IE_BS_2E_Angular_l)
{
    for (size_t L=0;L<=LMax;L++)
        this->Insert(new Orbital_IBS<K>(this,N,rmin,rmax,L));
    this->BuildCache(LMax);
}


template <size_t K> BasisSet<K>::BasisSet(size_t N, double rmin, double rmax, const ElectronConfiguration& ec)
: ::BSpline::BS_Common<K>(this, new IE_BS_2E_Angular_ml)
{
    const Atom_EC& aec=dynamic_cast<const Atom_EC&>(ec);
    size_t LMax=aec.GetLMax();
    for (size_t L=0;L<=LMax;L++)
    {
        auto mls=aec.GetBreadown(L);
        if (mls.ml_paired.size()>0)   
            ::BSpline::BS_Common<K>::Insert(new Orbital_IBS<K>(this,N,rmin,rmax,L,mls.ml_paired));            
        if (mls.ml_unpaired.size()>0)   
            ::BSpline::BS_Common<K>::Insert(new Orbital_IBS<K>(this,N,rmin,rmax,L,mls.ml_unpaired));            
        if (mls.ml_unoccupied.size()>0)   
            ::BSpline::BS_Common<K>::Insert(new Orbital_IBS<K>(this,N,rmin,rmax,L,mls.ml_unoccupied));            

    
    }           
    this->BuildCache(LMax);
}


#define INSTANCEk(k) template class BasisSet<k>;
#include "../../radial/BSpline/Instance.hpp"


}} //namespace

