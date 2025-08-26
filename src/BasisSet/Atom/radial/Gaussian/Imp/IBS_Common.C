// File: Atom/radial/Gaussian/IBS_Common.C  l/ml/kappa/mj independent part of Irrep Basis Set (IBS) for atom Gaussians.
module;
#include <iostream>
#include <iomanip>
#include <iosfwd>
#include <vector>
#include <cassert>
module qchem.BasisSet.Atom.Internal.radial.GaussianBS;
import qchem.BasisSet.Atom.Internal.radial.GaussianIntegrals;
import qchem.Symmetry;
import Common.IntPow;
import oml;



namespace Gaussian
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
    os << "Spherical Gaussian L=" << *GetSymmetry();
    itsEval->Write(os);
    return os;
}




} //namespace
