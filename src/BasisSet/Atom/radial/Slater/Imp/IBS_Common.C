// File: Atom/radial/Slater/IBS_Common.H  l/ml/kappa/mj independent part of Irrep Basis Set (IBS) for atom Slater functions.
module;
#include <iostream>
#include <cmath>
#include <vector>
#include <cassert>
module qchem.BasisSet.Atom.Internal.radial.SlaterBS;
import qchem.BasisSet.Atom.Internal.radial.Slater.Integrals;
import qchem.Symmetry;
import Common.IntPow;
import oml;

namespace Slater
{
//----------------------------------------------------------------
//
//  Common implementation for orbital and fit basis sets.
//
IrrepBasisSet::IrrepBasisSet(IBS_Evaluator* eval, Symmetry* sym)
    : IrrepBasisSet_Common<double>(sym)
    , itsEval(eval)
{};


std::ostream&  IrrepBasisSet::Write(std::ostream& os) const
{
    os << "Spherical Slater L=" << *GetSymmetry();
    itsEval->Write(os);
    return os;
}

} //namespace 