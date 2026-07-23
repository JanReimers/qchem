// File: BasisSet/Lattice_3D/BasisSet.C  Factory + container for 3D-periodic (plane-wave) basis sets.
//
// The public entry point for building a crystal basis set: hand it a Lattice_3D (the cell + its
// Brillouin-zone grid) and a cutoff, get back an abstract tBasisSet<dcmplx>.  The concrete container
// (PW_BasisSet) and the per-k PlaneWave_IBS list it owns are an implementation detail (Imp/BasisSet.C);
// callers are forced through the polymorphic BasisSet interface, exactly as for molecular bases.
module;
#include <memory>
export module qchem.BasisSet.Lattice_3D.BasisSet;
export import qchem.BasisSet;                          // Complex_BS (= tBasisSet<dcmplx>)
export import qchem.Lattice_3D;                        // Lattice_3D (the crystal structure + BZ grid)
export import qchem.BasisSet.Lattice_3D.PlaneWave_IBS; // PlaneWave_IBS + LocalPotential/SeparablePotential
export import qchem.BasisSet.Lattice_3D.GPW_IBS;       // GPW_IBS + CellImages (the GPWFactory mode argument)
import qchem.BasisSet.Internal.BasisSetImp;            // BasisSetImp<dcmplx> (the PW_BasisSet base; NOT re-exported)
import qchem.Types;                                    // dcmplx

export namespace qchem::BasisSet::Lattice_3D
{

//! \brief Which 3D-periodic basis to build.  PW = plane waves (lineage A); APW/LAPW (lineage B) follow.
enum class Type { PW };

//! \brief Build a 3D-periodic basis set for \a lat at energy cutoff \a Ecut.  Returns an abstract
//! tBasisSet<dcmplx> on the heap (caller owns), so callers must use the polymorphic interface.
//! \note The pseudopotential model is NOT configured here: it lives on the Vpseudo Hamiltonian term
//!   (the pseudo-wall), which calls the basis's MakeLocalPotential/MakeSeparablePotential assembly.
//! \note Single-k for now: the returned basis holds ONE Bloch block at \f$\Gamma\f$.  Phase 2
//!   generalises this to the full BZ k-list (intended to be the only k-loop in the framework).
Complex_BS* Factory(Type type, const ::qchem::Lattice_3D& lat, double Ecut);

//! \brief Build a GPW basis (periodic Gaussians on the lattice) over \a lat from a molecular Gaussian basis
//! \a mol built on the cell's atoms, at density-collocation cutoff \a densityEcut.
//! Returns an abstract tBasisSet<dcmplx> (caller owns).  Unlike the PW \c Factory this needs a Gaussian
//! orbital basis (GPW = Gaussian orbitals); the pseudopotential still lives on the Hamiltonian term, which
//! reaches GPW's real-space \c Integrals_Pseudo<dcmplx> assembly -- so the same \c Ham_PW_DFT drives it.
//! THERE IS NO CUT (doc/GPWPlan.md pin): every lattice sum is an eps-converged series enumerated inside the
//! molecular seam -- no radius parameter exists on this surface.
//! \param kShift  fractional Monkhorst-Pack offset of the k-mesh (\f$0\f$ = Γ-centred; \f$½\f$ = the classic MP
//!                offset, i.e. CP2K's default for even grids -- \f$k=\pm¼\f$ at \f$N=2\f$).
//! \param densityEcut  \f$<0\f$ = AUTOMATIC density grid \a cutoffFactor\f$\cdot\alpha_{\max}\f$ (recommended);
//!        \f$=0\f$ = 1E-only; \f$>0\f$ = explicit Hartree cutoff (\c cerr warning if under-resolved).
//! \param images  the lattice-image MODE (\c CellImages::Periodic default; \c HomeCellOnly = the box gates).
//! \param cutoffFactor  \f$C\ge4\f$ in the density-grid floor \f$C\cdot\alpha_{\max}\f$ (default 4).
Complex_BS* GPWFactory(const ::qchem::Lattice_3D& lat, std::shared_ptr<const BasisSet::Real_BS> mol,
                       double densityEcut, rvec3_t kShift={0,0,0},
                       CellImages images=CellImages::Periodic, double cutoffFactor=8.0);

} //namespace

// NOT exported: the concrete containers are an implementation detail (callers use Factory/GPWFactory's
// abstract Complex_BS).  Named here rather than anonymous in Imp/ so they are first-class types -- the home
// for the shared density-grid PeriodicGridEvaluator, GPW_BasisSet sitting beside PW_BasisSet.  Ctors in Imp/BasisSet.C.
namespace qchem::BasisSet::Lattice_3D
{

//! A tBasisSet<dcmplx> holding the plane-wave Bloch block(s); owns the IBS list (deleted with the basis).
class PW_BasisSet : public BasisSet::BasisSetImp<dcmplx>
{
public:
    //! Build one PlaneWave_IBS per Brillouin-zone k-point of \a lat at cutoff \a Ecut (a single Gamma block
    //! for an N=(1,1,1) mesh).  The basis ctor is the sole place that enumerates k -- the framework's
    //! per-irrep loop then IS the BZ sum \f$\sum_k w_k\f$.
    PW_BasisSet(const ::qchem::Lattice_3D& lat, double Ecut);
};

//! A tBasisSet<dcmplx> holding ONE \f$\Gamma\f$ GPW block (periodic Gaussians), built from a molecular
//! Gaussian basis over the cell's atoms + a density-collocation cutoff.  The GPW sibling of PW_BasisSet.
class GPW_BasisSet : public BasisSet::BasisSetImp<dcmplx>
{
public:
    GPW_BasisSet(const ::qchem::Lattice_3D& lat, std::shared_ptr<const BasisSet::Real_BS> mol,
                 double densityEcut, rvec3_t kShift={0,0,0},
                 CellImages images=CellImages::Periodic, double cutoffFactor=8.0);
};

} //namespace
