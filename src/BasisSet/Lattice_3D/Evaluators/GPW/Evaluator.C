// File: BasisSet/Lattice_3D/Evaluators/GPW/Evaluator.C  Gaussian-And-Plane-Waves orbital evaluator.
//
// The GPW sibling of PW_Evaluator (BasisSet/Lattice_3D/Evaluators/PW): a periodic ORBITAL evaluator whose
// orbitals are GAUSSIANS (compact, good for core/valence) instead of plane waves, standing on a lattice.  It
// satisfies the SAME isPW_1E_Evaluator concept, so the SAME EPW_Orbital1E_IBS<E> mixin builds the concrete
// GPW_IBS -- "GPW is a new evaluator, not a new IBS" (doc/MolecularPP_HarmonizationRound2.md section 2.5).
//
// The genuinely-new numerics of GPW's 1E tier are the PERIODIC Gaussian two-centre integrals -- the Bloch
// lattice sums  M_ij = Sum_R e^{ik.R} <chi_i | O | chi_j(.-R)>.  Rather than re-derive them, this evaluator
// DELEGATES them to the molecular Gaussian basis, which already owns the analytic McMurchie-Davidson kernels
// and GaussianRF::AtCenter (place a radial at an image): it holds a molecular orbital basis built over the
// cell's atoms and asks it the ONE high-level question Molecule::LatticeSum1E poses -- "sum these two-centre
// integrals over these lattice translations."  The Gaussian internals stay encapsulated on the molecular
// side; this evaluator only supplies the lattice geometry (the translation set from the cell) and complexifies
// the result.  No exponents or radials ever cross into qcLattice_BS.
//
// SCOPE (first increment): the GAMMA point (k=0).  The lattice sums are real (rsmat_t from the molecular side,
// widened to chmat_t here); general-k Bloch phases, the DFT tier (isPW_DFT_Evaluator: Hartree/XC by
// collocation) and a rigorous periodic nuclear attraction (Ewald) are later increments.
module;
#include <cstddef>
#include <functional>
#include <memory>
#include <string>
#include <vector>
export module qchem.BasisSet.Lattice_3D.Evaluators.GPW;
import qchem.BasisSet.Lattice_3D.Evaluators.PW;   // isPW_1E/DFT_Evaluator + PW_Grid_Evaluator (GPW's density grid)
import qchem.BasisSet.Molecule.LatticeSum1E;       // Molecule::LatticeSum1E (the periodic-1E capability we call)
import qchem.BasisSet;                             // Real_BS (the molecular Gaussian basis we own)
import qchem.BasisSet.Orbital_1E_IBS;              // Real_OIBS (its orbital block: op()/Gradient/size)
export import qchem.BasisSet.Internal.GMap;        // G_ERI3 (the DFT 3-centre tensor, now with GPW weights)
import qchem.UnitCell;                             // UnitCell (the direct lattice; CellsInSphere/ToCartesian)
import qchem.Structure;                            // Structure (the nuclear-attraction centres)
import qchem.Types;                                // rvec3_t, cvec_t, cvec3vec_t, chmat_t, rmat_t

export namespace qchem::BasisSet::Lattice_3D
{

//! \brief The GPW orbital evaluator: Gaussian orbitals on a lattice at a k-point (Gamma this increment).  Owns
//! the molecular Gaussian basis built over the cell's atoms and reaches its periodic-1E capability
//! (\c Molecule::LatticeSum1E) by an abstract cross-cast; supplies the lattice translation set from the cell.
//! Satisfies \c isPW_1E_Evaluator, so the \c EPW_Orbital1E_IBS mixin drives the concrete \c GPW_IBS.
class GPW_Evaluator
{
public:
    //! \param mol   the molecular Gaussian orbital basis built over \a cell's atoms (kept alive here); its one
    //!              orbital block must realise \c Molecule::LatticeSum1E (a PG Gaussian basis does).
    //! \param cell  the direct lattice (source of the real-space translation set \f$\{R\}\f$ + the reciprocal cell).
    //! \param densityEcut  plane-wave cutoff (Hartree) for the DENSITY/collocation grid: the FFT grid \f$\{r\}\f$
    //!              and the fit set \f$\{G\}\f$ GPW collocates its density onto (the DFT tier; must resolve the
    //!              Gaussian-product density -- CP2K's density cutoff).  \f$\le 0\f$ disables the DFT tier.
    //! \param kFrac fractional crystal momentum (this increment: \f$\Gamma\f$, asserted \f$\approx 0\f$).
    //! \param Rcut  radius (a.u.) of the lattice-translation sphere; \f$\le 0\f$ means the home cell only
    //!              (\f$R=0\f$: reproduces the finite molecule exactly).  A generous \a Rcut folds in images.
    GPW_Evaluator(std::shared_ptr<const BasisSet::Real_BS> mol, const UnitCell& cell,
                  double densityEcut = 0.0, const rvec3_t& kFrac = rvec3_t(0,0,0), double Rcut = 0.0);
    virtual ~GPW_Evaluator() = default;   // polymorphic: reached by the EPW_* mixin's Cast() cross-cast

    // --- isPW_1E_Evaluator surface (exact signatures the concept demands) ---
    size_t     size()                          const {return itsN;}
    cvec_t     Eval(const rvec3_t& r)          const;   //!< \f$\sum_R \chi_i(r-R)\f$ per function (Gamma: real)
    cvec3vec_t EvalGradient(const rvec3_t& r)  const;   //!< \f$\sum_R \nabla\chi_i(r-R)\f$
    chmat_t    OverlapMatrix()                 const;   //!< \f$\sum_R\langle i|j(\cdot-R)\rangle\f$
    chmat_t    KineticMatrix()                 const;   //!< \f$\sum_R\langle i|-\nabla^2|j(\cdot-R)\rangle\f$ (no 1/2)
    chmat_t    NuclearMatrix(const Structure* cl) const;//!< \f$\sum_R\langle i|\sum_c -Z_c/|r-R_c||j(\cdot-R)\rangle\f$

    // --- isPW_DFT_Evaluator surface: the DFT tier by COLLOCATION on the density grid (the genuinely-new GPW
    //     primitive).  All three route through the held PW_Grid_Evaluator (its {r}/{G}/FFT/Poisson). ---
    //! \brief Coulomb 3-centre tensor.  The per-column WEIGHT is the OVERLAP-metric (single-\f$r\f$) Fourier
    //! coefficient of the orbital product, \f$W_c(i,j)=\tfrac1\Omega\int\chi_i(r)\chi_j(r)\,e^{-iG_c\cdot r}\,d^3r\f$
    //! (\c ONE integration variable -- it is the collocated product FT, i.e. the density-side of the two-electron
    //! integral).  The SECOND electron coordinate and the \f$1/r_{12}\f$ live in the per-column diagonal Poisson
    //! kernel \f$4\pi/|G_c|^2\f$ (the \f$r_2\f$ integral \f$\int e^{iG_c\cdot r_2}/r_{12}\f$ in reciprocal space).
    //! So the full two-electron Coulomb is \c weight\f$\times\f$\c kernel, factorised through G-space -- never a
    //! 2-\f$r\f$ integral.  Density contracts \f$D\f$ against this via \c ContractG_ERI3.
    G_ERI3  Repulsion3CTensor() const;
    //! \brief Overlap 3-centre tensor: the same single-\f$r\f$ weight \f$W_c(i,j)\f$, EMPTY kernel -- the density's
    //! Fourier coefficient \f$\tilde\rho(G_c)=\sum_{ij}D_{ij}W_c(i,j)\f$ (no Poisson).
    G_ERI3  Overlap3CTensor() const;
    //! \brief The potential->KS-matrix bridge (collocation's adjoint): \f$\langle\chi_i|V|\chi_j\rangle=\int\chi_i
    //! V\chi_j\f$ with \f$V(r)\f$ the inverse-FFT of \a Vtilde over the density grid -- grid-integrate, not the PW
    //! Fourier lookup.  Satisfies \c isPW_DFT_Evaluator; forwarded by \c EPW_Orbital_DFT_IBS to \c MakeOverlap.
    chmat_t OverlapMatrix(const std::function<dcmplx(const ivec3_t&)>& Vtilde) const;

    //! The density/collocation grid engine (the fit basis is built over it, so \f$\tilde\rho\f$'s \f$\{G\}\f$ matches).
    const PW_Grid_Evaluator& DensityGrid() const {return *itsGrid;}

    //! Cache-key fragment: the molecular basis's ID + \f$k\f$ + translation count + the density-grid cutoff
    //! (the collocation tensor depends on the grid, so the framework cache key must pin it).
    std::string IDFragment() const;

private:
    std::shared_ptr<const BasisSet::Real_BS> itsMol;   //!< owns the molecular Gaussian basis (lifetime)
    const BasisSet::Real_OIBS*          itsOrb = nullptr; //!< its single orbital block (op()/Gradient/size)
    const Molecule::LatticeSum1E*       itsLat = nullptr; //!< the same block's periodic-1E capability (cross-cast)
    std::vector<rvec3_t>                itsR;             //!< Cartesian lattice translations (always incl. origin)
    rvec3_t                             itsk;             //!< fractional crystal momentum (Gamma this increment)
    size_t                              itsN   = 0;       //!< number of Gaussian orbitals
    std::shared_ptr<const PW_Grid_Evaluator> itsGrid;     //!< the density/collocation grid (null if DFT tier off)
    // NO hand-rolled tensor cache: the collocation tensor is a stateless build; the FRAMEWORK caches it
    // (BasisSet::Band_FT_IBS::Repulsion3C/Overlap3C via theCache<dcmplx>(), keyed by BasisSetID -- see IDFragment).
    rmat_t  PhiOnGrid()   const;   //!< \f$\chi_i(r_g)\f$ on the density grid (Npts x n; Gamma-real) -- computed on demand
    G_ERI3  BuildWeights() const;  //!< the collocation weight tensor (columns=grid {G}, weights, NO kernel)
};

static_assert(isPW_1E_Evaluator <GPW_Evaluator>, "GPW_Evaluator must satisfy isPW_1E_Evaluator");
static_assert(isPW_DFT_Evaluator<GPW_Evaluator>, "GPW_Evaluator must satisfy isPW_DFT_Evaluator");

} //namespace
