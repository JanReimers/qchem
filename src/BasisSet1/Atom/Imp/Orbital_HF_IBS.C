// File: BasisSet1/Atom/Imp/Orbital_HF_IBS.C 4 nested loops for HF integrals.
module;
// #include <iosfwd>
// #include <memory>
// #include <cassert>
// #include "forward.H"

module qchem.BasisSet.Atom.IBS;
// import qchem.BasisSet.Internal.IrrepBasisSetImp;
// import qchem.BasisSet.Internal.Orbital_DHF_IBS;
// import qchem.BasisSet.IrrepBasisSet;
// import qchem.BasisSet.Orbital_1E_IBS;
// import qchem.BasisSet.Orbital_DFT_IBS;
// import qchem.BasisSet.Orbital_HF_IBS;
import qchem.BasisSet.Atom.Evaluators.BS;
// import qchem.Symmetry.Yl;


template <isHF_Evaluator E> ERI4 Orbital_HF1_IBS::MakeExchange(const BasisSet::Orbital_HF_IBS<double>& _c) const 
{
    auto& a=dynamic_cast<const E&>(*this);
    auto& c=dynamic_cast<const E&>(_c);
    return itsEvaluator->Exchange(&a,&c);
}
