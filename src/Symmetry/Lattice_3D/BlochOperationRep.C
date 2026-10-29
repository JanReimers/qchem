// File: Symmetry/Lattice_3D/BlochOperationRep.C  A crystal op {W|tau} on a BLOCH orbital basis (doc/KAxisPlan.md P4).
//
// The lattice sibling of Molecule::BuildOperationRep: the same shell permutation x per-shell angular rep x per-component
// normalisation, plus what a crystal adds -- the image centre wraps back into the stored cell by a lattice vector t, and
// the Bloch sum turns that wrap into a PHASE.  With the cell-convention Bloch functions every lattice block in this code
// uses (GPW_Evaluator::Eval), \f$\phi_{k,a}(r)=\sum_R e^{ik\cdot R}\chi_a(r-A-R)\f$, which are periodic in k, an op
// \f$g\f$ maps block k onto block \f$k'=W_c k\f$ (any stored block with \f$k'\f$ mod G):
//
//     g phi_{k,Aa} = e^{-2 pi i k'.t} Sum_b (N_a/N_b) D_ba(W_c) phi_{k',A'b},     W A_f + tau = A'_f + t,
//
// with \f$k'_f=W^{-T}k_f\f$ fractional and \f$W_c=AWA^{-1}\f$.  Pure math on the shell list and the cell: no basis, no
// integrals.  The spin action sigma is NOT applied here -- it swaps the spin BLOCKS, which the caller owns.
// ⚠ That split ASSUMES COLLINEAR order: a spin op is then only a relabel of the up/down blocks.  Non-collinear support needs the
// SU(2) rep acting on a 2n spinor block (M (x) U(g), with time reversal antiunitary), and the callers that swap blocks change.
module;
#include <vector>
export module qchem.Symmetry.Lattice_3D.BlochOperationRep;
export import qchem.Symmetry.Molecule.OperationRep;   // AoShell (centre, angular rep, norms, offset)
export import qchem.Symmetry.SymOp;                   // SymOp {W|tau,sigma}, fractional direct frame
export import qchem.Types;                            // mat_t<dcmplx>, rvec3_t, rmat3d_t

export namespace qchem::Symmetry::Lattice_3D
{

//! \brief The fractional crystal momentum an op sends \a kFrac to: \f$k'=W^{-T}k\f$ (not reduced mod G).
rvec3_t ImageOfK(const SymOp& g, const rvec3_t& kFrac);

//! \brief The matrix \f$M\f$ carrying the Bloch orbitals at \a kFrac to those at \c ImageOfK(g,kFrac) under \a g:
//! \f$\hat g\,\phi_{k,\mu}=\sum_\nu M(\nu,\mu)\,\phi_{k',\nu}\f$, so orbital coefficients map as \f$C_{k'}=M C_k\f$.
//! \a shells: the block's shells (Cartesian centres, any angular convention); \a A: the cell matrix (columns = lattice
//! vectors); \a g in the fractional direct frame.  Unitary in the overlap metric, \f$M^\dagger S_{k'}M=S_k\f$, when \a g is
//! a symmetry of the structure AND the basis.  THROWS (std::runtime_error) when a shell has no image (\a g is not).
mat_t<dcmplx> BuildBlochOperationRep(const std::vector<Molecule::AoShell>& shells, const rmat3d_t& A,
                                     const SymOp& g, const rvec3_t& kFrac, double tol=1e-6);

} //namespace
