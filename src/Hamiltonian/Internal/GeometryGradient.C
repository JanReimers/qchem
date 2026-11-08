// File: Hamiltonian/Internal/GeometryGradient.C  The TERM face of the analytic force (doc/GeometryOptimizationPlan.md §13, §13.2; lane B
// 2026-11-08).  A term that can give its own geometry gradient implements it; tHamiltonianImp reaches it by the sanctioned abstract ->
// abstract cross-cast (as MakeResponseKernel reaches tResponse_HT), and a term WITHOUT it is a MissingGradients() entry: never a silent zero.
module;
#include <string>
export module qchem.Hamiltonian.Internal.GeometryGradient;
export import qchem.Hamiltonian;   // tbs_t, tDM_CD and the Types vocabulary
export import qchem.Structure;     // the atoms a force is taken in

export namespace qchem::Hamiltonian
{

using ChargeDensity::tStatic_PC;   // the per-block pullback client tDM_CD::DM_Pullback walks (the force twin of tStatic_CC)

template <class T> class tGeometryGradient_HT
{
public:
    virtual ~tGeometryGradient_HT() = default;
    //! EMPTY when this term can give its force on every block of \a bs for \a st; otherwise why not, as "<term>: <what is missing>" (the
    //! basis face that refuses, or the route).  CLIENT: the Hamiltonian's MissingGradients(), asked at setup before an optimization.
    virtual std::string MissingGradient(const tbs_t<T>* bs, const Structure* st) const=0;
    //! F[A] -= ∂E_term/∂R_A at the FIXED density \a D, per atom of \a st: the Hellmann-Feynman part with the basis motion included (the
    //! Pulay part, Tr(W ∂S), is the Hamiltonian's).  THROWS when a face it needs refuses: MissingGradient said so first.
    virtual void AddForces(rvec3vec_t& F, const tDM_CD<T>* D, const Structure* st) const=0;
};

} // namespace qchem::Hamiltonian
