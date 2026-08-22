// File: Hamiltonian/Types.C  Define types used throughout the Hamiltonian module.
module;

export module qchem.Hamiltonian.Types;

export import qchem.BasisSet;
export import qchem.BasisSet.Orbital_HF_IBS;
export import qchem.BasisSet.Orbital_DFT_IBS;
export import qchem.BasisSet.Orbital_DHF_IBS;
export import qchem.BasisSet.Orbital_1E_IBS;
export import qchem.BasisSet.Fit_IBS;
import qchem.Types;   // dcmplx (for cobs_t)


export namespace qchem::Hamiltonian
{
    template <class T> using tbs_t=BasisSet::tBasisSet<T>;         // whole (composite) basis: Iterate<tobs_t>() yields the per-irrep bases
    using rbs_t    =tbs_t<double>;
    using cbs_t    =tbs_t<dcmplx>;   // whole (composite) plane-wave basis (Ham_PW_DFT's fit-basis factory source)
    using fbs_t   =BasisSet::Fit_IBS;
    template <class T> using tobs_t=BasisSet::Orbital_1E_IBS<T>;  // T-parametric orbital basis
    // r* = <double>, c* = <dcmplx> (mirrors rsmat_t/chmat_t).
    using robs_t  =tobs_t<double>;  using cobs_t=tobs_t<dcmplx>;
    using rohfbs_t =BasisSet::Orbital_HF_IBS<double>;
    using odftbs_t=BasisSet::Orbital_DFT_IBS<double>;
    using orkbbs_t=BasisSet::Orbital_RKB_IBS<double>;

//! \brief WHICH FIT BASIS represents \f$v_{xc}\f$ (doc/SymmetryUpgradePlan.md §6a, the fit/grid
//! SEPARATION, user 2026-08-01): the fit-basis choice and the real-space grid choice
//! (\c qcMesh::MeshParams) are ORTHOGONAL user knobs.
//!  - \c PlaneWave: expand \f$v_{xc}\f$ on the \f$\{Q_j\}\f$ ball (band-limited; the projection
//!    quadrature is the FFT on the uniform raster) -- the PAIR/collocation \c XC_Quadrature.
//!  - \c Delta: the delta-function "fit" -- coefficients ARE the grid-point values, H by direct
//!    quadrature -- the SINGLES (Φ-table) \c XC_Quadrature, on ANY real-space grid (Becke or uniform).
//!  - \c Auto: picks Delta whenever the plane-wave fit cannot do the job -- on a Becke grid (no G-space
//!    raster) and on a POLARIZED run (the pair route is not spin-native) -- else PlaneWave.
//! (PlaneWave fit ON a Becke grid = I3: the projection sum is trivial, but the one-functional E/H
//! derivative pairing must be designed first -- asserted out until then.)
//!
//! PUBLIC (moved out of qchem.Hamiltonian.Internal.Hamiltonians 2026-08-08): it is a user-facing policy
//! knob, so a FACADE must be able to name it without importing an internal across a library boundary.
//!
//! 2026-08-22: the ENUM itself moved DOWN to \c qchem.BasisSet.Fit_IBS, beside the fit-basis types it
//! selects between, because the factory that must read it -- \c CreateVxcFitBasisSet -- lives in
//! qcBasisSet and cannot see qcHamiltonian.  This alias keeps the user-facing spelling
//! \c Hamiltonian::VxcFit::Delta exactly as it was; \c Ham_PW_DFT resolves \c Auto (the one decision that
//! needs to know the run is polarized) and passes the answer to the factory.
using VxcFit = BasisSet::VxcFit;
}
