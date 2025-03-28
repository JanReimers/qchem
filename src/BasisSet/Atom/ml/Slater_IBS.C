// File: Atom/ml/Slater_IBS.C  Slater Irrep Basis Set (IBS) with orbital angular momentum l,m.

#include "Imp/BasisSet/Atom/ml/Slater_IBS.H"
#include "Imp/BasisSet/Atom/ml/Slater_BF.H"
#include "Imp/BasisSet/Atom/radial/Slater/Integrals.H"
#include "Imp/Symmetry/YlmQN.H"
#include <iostream>
#include <cassert>

namespace Atom_ml
{
namespace Slater
{



::Fit_IBS* Orbital_IBS::CreateCDFitBasisSet(const ::BasisSet*,const Cluster*) const
{
    // return new IrrepBasisSet(itsLAParams,GetDataBase(),0,es*2,0,0);
    assert(false);
    return 0;
}

::Fit_IBS* Orbital_IBS::CreateVxcFitBasisSet(const ::BasisSet*,const Cluster*) const
{
    // return new IrrepBasisSet(itsLAParams,GetDataBase(),0,es*2.0/3.0,0,0);
    assert(false);
    return 0;
}

// std::ostream&  IrrepBasisSet::Write(std::ostream& os) const
// {
//     if (Pretty())
//     {
//         os << "Slater functions l,m=" << GetQuantumNumber()
//         << " with " << GetNumFunctions() << " basis functions, alpha={";
//         for (auto b:*this) os << *b;
//         os << "}" << std::endl;
//     }
//     return os;
// }



::IrrepBasisSet* Orbital_IBS::Clone(const RVec3&) const
{
    std::cerr << "Why are you relocating a Slater atomic basis set?!" << std::endl;
    return 0;
}


}} //namespace
