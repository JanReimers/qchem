// File: Atom/l/Gaussian_IBS.H  Gaussian Irrep Basis Set (IBS) with orbital angular momentum l.

#include "Imp/BasisSet/Atom/l/Gaussian_IBS.H"
#include "Imp/BasisSet/Atom/radial/Gaussian/BasisFunction.H"
#include "Imp/BasisSet/Atom/l/GaussianIE.H"
#include "Imp/BasisSet/Atom/radial/Gaussian/Integrals.H"
#include "Imp/Symmetry/YlQN.H"
#include <BasisSet.H>
#include <iostream>
#include <cassert>

namespace Atoml
{
namespace Gaussian
{
  

//----------------------------------------------------------------
//
// Orbital SG basis set.
//

::Fit_IBS* Orbital_IBS::CreateCDFitBasisSet(const ::BasisSet* bs,const Cluster*) const
{
    auto db=dynamic_cast<const DB_cache<double>*>(bs);
    return new ::Gaussian::Fit_IBS(db,es*2,0);
}

::Fit_IBS* Orbital_IBS::CreateVxcFitBasisSet(const ::BasisSet* bs,const Cluster*) const
{
    auto db=dynamic_cast<const DB_cache<double>*>(bs);
    return new ::Gaussian::Fit_IBS(db,es*2.0/3.0,0);    
}

::IrrepBasisSet* Orbital_IBS::Clone(const RVec3&) const
{
    std::cerr << "Why are you relocating a spherical Gaussian basis set?!" << std::endl;
    return 0;
}



} //namespace
} //namespace
