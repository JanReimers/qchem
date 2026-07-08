// File: BasisSet/Molecule/LatticeSum1E.C  Periodic (lattice-summed) 1-electron integrals of a molecular
// Gaussian basis -- the GPW seam.
//
// A molecular Gaussian orbital basis is a set of contracted Gaussians standing at the atoms.  Placed in a
// periodic cell, the Bloch orbital is the lattice sum  chi^k_i(r) = Sum_R e^{ik.R} chi_i(r-R), and its
// one-electron matrices are the corresponding lattice sums of the ordinary (finite) two-centre integrals:
//   S_ij = Sum_R e^{ik.R} <chi_i | chi_j(.-R)> ,  and likewise <p^2> and the nuclear attraction.
// The R=0 term IS the finite integral, so with R={0} these reproduce the molecule's own Overlap/Kinetic/
// Nuclear matrices exactly.
//
// This is the ONE high-level operation GPW needs from the molecular basis: it hands the basis a list of
// real-space lattice translations (Cartesian, MUST include {0,0,0}) and gets back the summed matrices.  The
// Gaussian machinery -- the primitives, the M&D kernels, GaussianRF::AtCenter that places a radial at each
// image -- stays ENCAPSULATED on the molecular side (the caller never sees a radial or an exponent); GPW
// stays purely lattice-side.  A concrete molecular Gaussian basis realises this face and GPW reaches it by an
// abstract->abstract cross-cast of the orbital block.
//
// ENGINE-NEUTRAL SEAM (the integral-engine switch point).  This face is deliberately engine-agnostic: it
// names no integral engine, so it can be realised by ANY molecular Gaussian basis.  Today the McMurchie-
// Davidson Cartesian basis (PG_Cart::Orbital_IBS, Engine::MnD) implements it -- its AtCenter + analytic 2C
// kernels make the per-image sum exact and trivial, ideal for the correctness increment.  The FASTER libCint
// engine (PG_LibCint) can realise the same face (per image, the cross-centre integrals between the home and
// R-shifted shells) -- a perf follow-up.  Because GPW consumes ONLY this abstract face (it takes any
// Molecule::Real_BS and cross-casts, throwing cleanly if unsupported), switching engines is then just the
// Engine argument to Molecule::Factory where the basis is built -- NO change to GPW itself.
//
// SCOPE (first increment): the GAMMA point (k=0), where every e^{ik.R}=1 and the sums are REAL -- hence the
// rsmat_t returns and no phase argument.  The general-k Bloch case (complex phases e^{ik.R} -> chmat_t) and a
// physically rigorous periodic nuclear attraction (Ewald, vs. this large-cell single-image limit) are later
// increments; this contract is expected to grow with them.
module;
#include <vector>
export module qchem.BasisSet.Molecule.LatticeSum1E;
import qchem.Structure;   // Structure (the nuclear-attraction centres)
import qchem.Types;       // rvec3_t
import qchem.Blaze;       // rsmat_t

export namespace qchem::BasisSet::Molecule
{

//! \brief The periodic (lattice-summed) 1-electron capability of a molecular Gaussian basis: the GAMMA-point
//! Bloch matrices \f$M_{ij}=\sum_R\langle\chi_i|\,\hat O\,|\chi_j(\cdot-R)\rangle\f$ over a caller-supplied
//! set of Cartesian lattice translations \a Rs (which MUST contain the origin).  With \c Rs={0} each returns
//! the finite molecule's own matrix; growing \a Rs folds in the periodic images.  Real at \f$\Gamma\f$.
//!
//! \c Make prefix / not self-caching: like the other 1E builders, these are the single-shot producers behind
//! the cached \c Integrals_{Overlap,Kinetic,Nuclear} accessors (\c Overlap()/\c Kinetic()/\c Nuclear() on the
//! GPW orbital IBS cache them in \c DB_Cache keyed on \c BasisSetID, so a static-across-SCF matrix is built
//! once and reused).  The \a Rs argument already marks these as the lattice (periodic) forms -- no name prefix.
class LatticeSum1E
{
public:
    virtual ~LatticeSum1E() = default;

    //! \f$S_{ij}=\sum_R\langle\chi_i|\chi_j(\cdot-R)\rangle\f$ (normalised).
    virtual rsmat_t MakeOverlap(const std::vector<rvec3_t>& Rs) const = 0;
    //! \f$\langle p^2\rangle_{ij}=\sum_R\langle\chi_i|-\nabla^2|\chi_j(\cdot-R)\rangle\f$ (NO 1/2 -- the
    //! Hamiltonian applies it, matching \c Integrals_Kinetic).
    virtual rsmat_t MakeKinetic(const std::vector<rvec3_t>& Rs) const = 0;
    //! \f$\sum_R\langle\chi_i|\sum_c -Z_c/|r-R_c|\,|\chi_j(\cdot-R)\rangle\f$ over the atoms of \a cl.
    virtual rsmat_t MakeNuclear(const std::vector<rvec3_t>& Rs, const Structure* cl) const = 0;
};

} //namespace
