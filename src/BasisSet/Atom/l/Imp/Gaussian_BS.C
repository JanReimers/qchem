// File: Atom/l/Gaussian_BS.H
module;
#include <cmath>
#include <vector>
#include <iostream>
module qchem.BasisSet.Atom.Internal.l.GaussianBS;
import qchem.BasisSet.Atom.Internal.radial.Gaussian.ExponentScaler; 
import qchem.Symmetry.AtomEC;

namespace Atoml
{
namespace Gaussian
{

BasisSet::BasisSet(size_t N, double emin, double emax, size_t LMax)
{
    ::Gaussian::ExponentScaler gs(N,emin,emax,LMax);
    for (size_t L=0;L<=LMax;L++)
        Insert(new Orbital_IBS(this,gs.Get_es(L),L)); 
}

} //namespace
} //namespace

namespace Atom_ml
{
namespace Gaussian
{
BasisSet::BasisSet(size_t N, double emin, double emax, const ElectronConfiguration& ec)
{
    const Atom_EC& aec=dynamic_cast<const Atom_EC&>(ec);
    ::Gaussian::ExponentScaler ss(N,emin,emax,aec.GetLMax());
    for (size_t L=0;L<=aec.GetLMax();L++)
    {
        auto mls=aec.GetBreadown(L);
        if (mls.ml_paired.size()>0)   
            Insert(new Atoml::Gaussian::Orbital_IBS(this,ss.Get_es(L),L,mls.ml_paired));            
        if (mls.ml_unpaired.size()>0)   
            Insert(new Atoml::Gaussian::Orbital_IBS(this,ss.Get_es(L),L,mls.ml_unpaired));            
        if (mls.ml_unoccupied.size()>0)   
            Insert(new Atoml::Gaussian::Orbital_IBS(this,ss.Get_es(L),L,mls.ml_unoccupied));            
    }
}

} //namespace
} //namespace
