// File: BasisSet/Band_FT_IBS.C  Abstract reciprocal-space (Fourier / plane-wave) DFT-assembly capability.
//
// The dual of Band_DFT_IBS: instead of integrating real-space ScalarFunctions on a mesh, this is the
// RECIPROCAL-SPACE lineage -- a plane-wave basis and a plane-wave density share G-space coefficients (a
// ΔG_Map, rho-tilde).  The questions here are the DENSITY-driven assembly: the density-fit-basis factory,
// the D-free {G} 3-centre gather (the density contracts D against it, off-basis), and the FFT-free Hartree
// matrix from rho-tilde.  (The FFT XC grid
// engine lives on the BasisSet::G_FieldEvaluator seam; the EXTERNAL pseudopotential assembly is a separate
// capability, Pseudopotential::Integrals_Pseudo in qcPseudopotential, so qcBasisSet names no pseudopotential type.)
// A term reaches this the sanctioned way: holding the abstract orbital basis and dynamic_cast-ing UP.
module;
#include <functional>
export module qchem.BasisSet.Band_FT_IBS;
export import qchem.BasisSet.Orbital_1E_IBS;
export import qchem.BasisSet.Internal.GMap;
import qchem.Types;       // hmat_t<dcmplx>
import qchem.BasisSet.Fit_IBS;   // cFIT_CD_ABS (the auxiliary density-fit basis this basis creates) + qcMesh::MeshParams

export namespace qchem::BasisSet
{

//! \brief A basis that assembles matrices in reciprocal (G-)space: the density-driven KS matrices from a
//! density's Fourier coefficients, and the external (pseudo)potential matrices from a per-species form
//! factor.  Only the complex (plane-wave) path realizes this; a density exposes its rho-tilde via the
//! matching ChargeDensity FourierDensity sub-interface, and the BZ sum is the per-block rho-tilde summed.
class Band_FT_IBS
    : public virtual Orbital_1E_IBS<dcmplx>
{
public:
    //! \brief Create THIS basis's auxiliary density-fit basis -- the plane-wave analog of
    //! Orbital_DFT_IBS::CreateCDFitBasisSet.  A Hartree/DFT term obtains its density fitter THROUGH this
    //! factory (never assuming orbital==fit): a distinct \c cFIT_CD_ABS over the tunable \f$\{G\}\f$ grid
    //! (= the orbital grid now; denser for the density later).  \a mp is accepted for interface symmetry
    //! with the molecular factory (a plane-wave fit basis is \f$E_{cut}\f$/grid-based, so it is ignored).
    //! Caller owns the result.
    virtual cFIT_CD_ABS* CreateCDFitBasisSet(const Structure* cl, const qcMesh::MeshParams& mp) const=0;

    //! \brief Create THIS basis's auxiliary potential (Vxc) fit basis -- the overlap-metric sibling of
    //! CreateCDFitBasisSet (the plane-wave analog of Orbital_DFT_IBS::CreateVxcFitBasisSet).  A distinct
    //! \c cFIT_SF_ABS over the tunable \f$\{G\}\f$ grid (same grid as CD today; the two diverge with the
    //! future denser-\f$\{G\}\f$ upgrade).  Caller owns the result.
    virtual cFIT_SF_ABS* CreateVxcFitBasisSet(const Structure* cl, const qcMesh::MeshParams& mp) const=0;

    //! \brief D-free COULOMB three-centre tensor: the delta support with the diagonal Poisson kernel
    //! \f$4\pi/|G_c|^2\f$ filled.  A density contracts \f$D\f$ against it (\c ContractG_ERI3) to get \f$V_H\f$
    //! directly.  Mirrors \c Orbital_DFT_IBS::Repulsion3C (Coulomb metric); \a c is the CD fit basis.  Cached
    //! (the process-wide integrals cache, keyed by BasisSetID); the one-time build is \c MakeRepulsion3C.
    const G_ERI3& Repulsion3C(const cFIT_CD_ABS& c) const;
    //! \brief D-free OVERLAP three-centre tensor (delta support, EMPTY kernel).  Mirrors
    //! \c Orbital_DFT_IBS::Overlap3C (overlap metric); \a c is the Vxc fit basis.  Cached; build is \c MakeOverlap3C.
    const G_ERI3& Overlap3C(const cFIT_SF_ABS& c) const;

    // NB: the SAD seed's structure-factor density (MakeFourierDensity) is NOT here -- it is a grid-engine
    // operation (G_FieldEvaluator::MakeFourierDensity), which the seed reaches through its OWN fit basis, so
    // the seed never depends on the orbital basis.

    // NB: the Hartree matrix is NOT built here.  A density contracts D against Repulsion3C to a V_H (kernel
    // baked), and the term assembles <i|V_H|j> via the G_FieldEvaluator::MakePotential primitive -- so there is
    // no rho-tilde->matrix method on this face (the old Repulsion(ΔG_Map) fitter path is retired).

    // NB: the FFT XC route -- the real-space grid, the inverse/forward transforms, and the
    // <i|V|j>=Vtilde(m_i-m_j) assembly -- is NOT here.  It is the plane-wave grid engine, exposed by the
    // BasisSet::G_FieldEvaluator seam (implemented by PW_Evaluator, so both the orbital and the auxiliary fit
    // basis carry it); the Vxc term quadratures on its FIT basis's grid through that seam, not on this one.

protected:
    //! \brief One-time build of the Coulomb tensor (the concrete plane-wave basis's own \f$\{G\}\f$ delta
    //! support + \f$4\pi/|G_c|^2\f$ kernel).  Called once per (basis, fit-basis) key by the cached
    //! \c Repulsion3C above -- the plane-wave analogue of \c Orbital_DFT_IBS::MakeRepulsion3C.
    virtual G_ERI3 MakeRepulsion3C(const cFIT_CD_ABS& c) const=0;
    //! \brief One-time build of the overlap tensor (delta support, empty kernel).  Called once by the cached
    //! \c Overlap3C -- the plane-wave analogue of \c Orbital_DFT_IBS::MakeOverlap3C.
    virtual G_ERI3 MakeOverlap3C(const cFIT_SF_ABS& c) const=0;
};

} //namespace
