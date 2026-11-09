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
    //! THE CONTRACT WITH THE PULAY SUM: W is built from the eigen-system of the terms' MATRICES, so a term's force is taken with ITS MATRIX
    //! as its D-gradient.  Where a term's energy depends on D also through a piece its matrix leaves out, that piece's R-derivative is
    //! left out too, because the Pulay sum cancels it.  The one case so far is a charge-constrained fit, whose target N = Tr(DS) moves
    //! with the basis.  Its λ Tr(D ∂S/∂R) cancels the λ D that a matrix V + λS would add to W.  λS shifts every level by λ and moves no
    //! orbital, so the matrix leaves it out (Vee_DFT_Dunlap).  The per-term fixed-D finite difference of such a term differs from its
    //! force by exactly that piece, so its gate is the whole force against the re-converged energy.
    virtual void AddForces(rvec3vec_t& F, const tDM_CD<T>* D, const Structure* st) const=0;
};

} // namespace qchem::Hamiltonian
