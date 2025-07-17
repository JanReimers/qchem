// File Slater_m/BasisSet.H

#include <iostream>
#include <vector>
#include "ml/Slater_BS.H"
#include "ml/Slater_IBS.H"
#include "radial/Slater/ExponentScaler.H"
#include "Symmetry/Atom_EC.H"
#include "Common/stl_io.h"

using std::cout;
using std::endl;

namespace Atom_ml
{
namespace Slater
{

BasisSet::BasisSet(size_t N, double emin, double emax, const ElectronConfiguration& ec)
{
    const Atom_EC& aec=dynamic_cast<const Atom_EC&>(ec);
    ::Slater::ExponentScaler ss(N,emin,emax,aec.GetLMax());
    for (size_t L=0;L<=aec.GetLMax();L++)
    {
        auto mls=aec.GetBreadown(L);
        if (mls.ml_paired.size()>0)   
            Insert(new Orbital_IBS(this,ss.Get_es(L),L,mls.ml_paired));            
        if (mls.ml_unpaired.size()>0)   
            Insert(new Orbital_IBS(this,ss.Get_es(L),L,mls.ml_unpaired));            
        if (mls.ml_unoccupied.size()>0)   
            Insert(new Orbital_IBS(this,ss.Get_es(L),L,mls.ml_unoccupied));            

    
    }
    
}

}} //namespace
