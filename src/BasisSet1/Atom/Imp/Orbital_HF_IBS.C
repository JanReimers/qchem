// File: BasisSet1/Atom/Imp/Orbital_HF_IBS.C 4 nested loops for HF integrals.
module;
// #include <iosfwd>
// #include <memory>
// #include <cassert>
// #include "forward.H"

module qchem.BasisSet1.Atom.IBS;
// import qchem.BasisSet1.Internal.IrrepBasisSetImp;
// import qchem.BasisSet1.Internal.Orbital_DHF_IBS;
// import qchem.BasisSet1.IrrepBasisSet;
// import qchem.BasisSet1.Orbital_1E_IBS;
// import qchem.BasisSet1.Orbital_DFT_IBS;
// import qchem.BasisSet1.Orbital_HF_IBS;
import qchem.BasisSet1.Atom.Evaluators.BS;
// import qchem.Symmetry.Yl;


template <isHF_Evaluator E> ERI4 Orbital_HF1_IBS::MakeExchange(const BasisSet1::Orbital_HF_IBS<double>& _c) const 
{
    auto& a=dynamic_cast<const E&>(*this);
    auto& c=dynamic_cast<const E&>(_c);
    return itsEvaluator->Exchange(&a,&c);
}
