// File: Atom/ml/BSpline_BS.C BSpline Basis Set for atoms, no m degeneracy.
module;
#include <iostream>
#include <cassert>

module qchem.BasisSet.Atom.ml.BSplineBS;
import qchem.Basisset.Atom.radial.BSpline.Rk;
import qchem.BasisSet.Imp.Cache4;
import qchem.Symmetry.AtomEC;
namespace Atom_ml
{
namespace BSpline
{


template <size_t K> BasisSet<K>::BasisSet(size_t N, double rmin, double rmax, const ElectronConfiguration& ec)
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
#include "../../Instance.hpp"
}} //namespace
