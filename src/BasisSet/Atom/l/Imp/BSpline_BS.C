// File: Atom/l/BSpline_BS.C BSpline Basis Set for atoms.
module;
#include <iostream>
#include <cassert>

module qchem.BasisSet.Atom.Internal.l.BSplineBS;
import qchem.Basisset.Atom.radial.BSpline.Rk;
import qchem.BasisSet.Internal.Cache4;
import qchem.BasisSet.Atom.Internal.Angular;

import qchem.Symmetry.AtomEC;

namespace Atoml
{
namespace BSpline
{


template <size_t K> BasisSet<K>::BasisSet(size_t N, double rmin, double rmax, size_t LMax)
: ::BSpline::BS_Common<K>(new IE_BS_2E_Angular_l)
{
    for (size_t L=0;L<=LMax;L++)
        this->Insert(new Orbital_IBS<K>(this,this,N,rmin,rmax,L));
    this->BuildCache(LMax);
}


template <size_t K> BasisSet<K>::BasisSet(size_t N, double rmin, double rmax, const ElectronConfiguration& ec)
: ::BSpline::BS_Common<K>(new IE_BS_2E_Angular_ml)
{
    const Atom_EC& aec=dynamic_cast<const Atom_EC&>(ec);
    size_t LMax=aec.GetLMax();
    for (size_t L=0;L<=LMax;L++)
    {
        auto mls=aec.GetBreadown(L);
        if (mls.ml_paired.size()>0)   
            ::BSpline::BS_Common<K>::Insert(new Orbital_IBS<K>(this,this,N,rmin,rmax,L,mls.ml_paired));            
        if (mls.ml_unpaired.size()>0)   
            ::BSpline::BS_Common<K>::Insert(new Orbital_IBS<K>(this,this,N,rmin,rmax,L,mls.ml_unpaired));            
        if (mls.ml_unoccupied.size()>0)   
            ::BSpline::BS_Common<K>::Insert(new Orbital_IBS<K>(this,this,N,rmin,rmax,L,mls.ml_unoccupied));            

    
    }           
    this->BuildCache(LMax);
}

template <size_t K> Vector<double> BasisSet<K>::loop_4_direct(size_t id, size_t la, size_t lc)  const
{
    Vector<double> rk(1);
    const Cacheable* c=Cache4::loop_4(id);
    const ::BSpline::RkEngine<K>* cd = dynamic_cast<const ::BSpline::RkEngine<K>*>(c);
    rk(1)= cd->Coulomb_R0();
    return rk;
}

#define INSTANCEk(k) template class BasisSet<k>;
#include "../../radial/BSpline/Instance.hpp"


}} //namespace

