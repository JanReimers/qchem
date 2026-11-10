// File: Calculation/SolidState.C  A periodic SCF state ON DISK: the SOLID's half -- its fingerprint (cell, k-blocks, GPW knobs, +U),
// its block walk, and the rebuild of the density it carries (doc/OpenWork.md §2 "SCF checkpoint/restart", CK-1).  What every kind
// shares (the framing, the writer and reader, the comparison, the Löwdin carry) is \c qchem.SCFState since 2026-11-09
// (V-FACADE-STATE); this module re-exports it.
//
// WHY IT EXISTS.  Three needs, one mechanism: crash-resume (the Oct 6-20 unattended window), a LIBRARY of
// converged states per material (user 2026-09-27: "restart any material from a good state"), and the
// self-consistent-U outer loop picking up across PROCESSES where the U = 0 run finished.  In-process the
// facade already warm-starts every outer step from its own density; this is the across-process half.
//
// WHAT IS SAVED, per (k,σ) block -- our analogue of QE's `.save` (inspected 2026-09-27: per-k wfc, ρ(G),
// occup.txt, the input echo):
//   D     the block's AO density matrix, PHYSICAL occupations (the BZ weight is stored beside it, never in it)
//   C,ε,f ALL orbitals (occupied and virtual), for CK-2 (a WaveFunction read from disk: χ₀/ACBN0/gaps on a
//         stored state with no SCF).  CK-1 restores from D alone.
// plus the FINGERPRINT and the run's own summary (energy, converged, iterations, the last stage's recipe).
//
// WHAT IS DELIBERATELY NOT SAVED.  No mixer or accelerator history: a restart is a FRESH STAGE (the
// 2026-09-21 stale-Pulay lesson -- eight near-zero residuals extrapolate the first new step straight back
// onto the old density).  No +U occupation matrix: ours is a FUNCTION of D (occupations from D_out), so it
// is recomputed from the restored D, never a second source of truth.  (Frozen-occupation mode -- LR,
// polarons -- would make n independent state; it does not exist yet, and when it does it must be saved.)
//
// THE FINGERPRINT has two classes, and the split IS the design:
//   REFUSE  the saved D is a matrix over THESE functions at THESE k-points for THIS many electrons: the
//           cell, the atom kinds, the species/PP (set + q per element), the orbital basis's shells, the k-mesh
//           blocks, Nelec, the spin group.  Any difference and the state means nothing here -- a REFUSED
//           Outcome, never a silently wrong start.
//   WARM    a restart ACROSS these is exactly what the mechanism is for: +U (the SCU loop), the density
//           grids and XC mesh (grid continuation), the multiplicity, the orthogonaliser, kT and the recipe
//           (anneals), a block's real/complex working type, the FUNCTIONAL (LDA -> PBE), and the atom POSITIONS
//           with the basis centres that follow them (a geometry step; the max displacement is printed).
//           No difference at all = an EXACT RESUME.
//   (C4 step 5: positions, shell centres and the functional used to be REFUSE; a state written under that
//   schema is migrated on comparison, so it still loads.)
//
// The file is read in Python by h5py (the Viz plan's carrier): complex arrays arrive as complex128, and the
// layout is written out below the classes.
module;
#include <map>
#include <memory>
#include <string>
#include <vector>
export module qchem.SolidState;
export import qchem.SCFState;            // the shared half: RestartRefusal, StateFingerprint, SavedState, ReadState, LowdinCarry...
import qchem.Types;
import qchem.Structure;
import qchem.BasisSet;                   // Complex_BS (the Bloch blocks)
import qchem.WaveFunction;               // cWaveFunction (the orbitals a save reads)
import qchem.ChargeDensity;              // cDM_CD (the restored density)
import qchem.Hamiltonian.Factory;        // HubbardManifold, SpinGroup

export namespace qchem
{

//! The \c grids.raster value of a run on the ball's own rasters (every raster demand 0).  States written before 2026-11-04 say
//! "BallOnly" for the same thing; reading one maps it here, so they resume as before.
inline constexpr const char* BallRastersLabel="ballRasters";

//! The inputs of the solid's own keys that only the FACADE knows (its options); the rest is asked of the built basis.  (Nelec, the
//! multiplicity, the functional and kT are the shared \c Calculation::StateIdentity's since 2026-11-09.)
struct RunIdentity
{
    std::vector<std::pair<std::string,int>>      species;      //!< (element, valence)
    std::string                                  ppSet;        //!< the pseudopotential SET's name (the family; q per Z is in \c species)
    std::vector<Hamiltonian::HubbardManifold>    hubbard;      //!< the U the Hamiltonian CURRENTLY carries
    double                                       densityEcut = 0.0, cutoffFactor = 0.0;
    std::string                                  xcMesh;       //!< the RESOLVED XC quadrature, named
    int                                          ortho = 0;
    double                                       orthoTol = 0.0;
    //! The orbital basis FAMILY, "gaussian" or "planeWave" (REFUSE: D is a matrix over one family's functions).  A Gaussian
    //! basis is identified by its shells (AoShellSource); a plane-wave basis by \c planeWaveEcut -- with the cell and each
    //! block's k and size (both already REFUSE) that fixes its G-set.
    std::string                                  orbitalBasis = "gaussian";
    double                                       planeWaveEcut = 0.0;   //!< Ha, the plane-wave orbital cutoff (planeWave only)
    //! \name The numerics that move the converged ground state without changing what D is a matrix OVER (all WARM; added
    //! 2026-10-30, OpenWork §3 "A restart's fingerprint is incomplete" gap 2).  A state written before them compares as agreeing.
    //!@{
    std::string                                  vxcFit;       //!< the RESOLVED v_xc fit basis ("Delta" | "PlaneWave")
    double                                       beckeEps = 0.0;   //!< the periodic Becke competitor-series eps (Becke mesh only, else 0)
    std::string                                  raster;       //!< what the rasters were sampled on: \c BallRastersLabel, or the demands (Gaussian family only)
    double                                       ladderFactor = 0.0;   //!< the multigrid ladder gradation (Gaussian family only)
    //! The GPW tolerances by name (vlocEps, screenEps, relCutoff, mgridEcuts, ...; Gaussian family only), each recorded as
    //! \c gpw.<name>: the facade knows the typed struct, this module only the values.
    std::map<std::string, std::vector<double>>   gpwTolerances;
    //!@}
    //! The RESOLVED CP2K-deviation routes that move the converged state (\c Deviation::movesState), keyed by their deck knob
    //! (\c policy.<route>) -> on/off; all WARM, recorded one key per route so a banner reformat is not a difference (added
    //! 2026-10-30, OpenWork §3 "A restart's fingerprint does not record the run POLICY").  A state without them compares as agreeing.
    std::map<std::string, bool>                  policy;
};

//! The SOLID's own fingerprint keys, beside the shared ones \c Calculation::RecordState writes (structure, electrons, blocks, scf.kT):
//! the cell matrix \a A (columns = lattice vectors), the species and PP set, the orbital basis of the built Bloch basis \a bs (the
//! shells off the AoShellSource face, or a plane-wave cutoff), each block's k and BZ weight under spin group \a g, and the facade's
//! \a id (+U, grids, GPW tolerances, policy routes).
void AddSolidKeys(StateFingerprint& fp, const rmat3d_t& A, const BasisSet::Complex_BS& bs, SpinGroup g, const RunIdentity& id);

//! \brief The solid's blocks in the canonical order (spin irreps outer, basis blocks inner -- the order tCompositeWF builds its
//! children in), each with its k as the block extra "k" and its BZ weight; with their orbitals from \a wf when given.
std::vector<StateBlock> SolidStateBlocks(const BasisSet::Complex_BS& bs, SpinGroup g, const WaveFunction::cWaveFunction* wf=nullptr);

//! \brief The solid's schema history, applied to a SAVED fingerprint before it is compared with this run's \a now: a state written
//! before C4 step 5 filed the functional, the atom positions and the shell centres under REFUSE; keys added since (pp.set,
//! basis.family, the numerics, the rung gate, the raster demands, the policy routes, the block irrep labels) take this run's value.
StateFingerprint MigrateSolidFingerprint(StateFingerprint saved, const StateFingerprint& now);
//! \c CompareFingerprints under \c MigrateSolidFingerprint: every solid comparison goes through here, so an old state keeps loading.
Outcome<std::vector<std::string>, RestartRefusal> CompareSolidFingerprints(const StateFingerprint& saved, const StateFingerprint& now);

//! The solid's in-memory snapshot (\c SnapshotState over \c SolidStateBlocks): the carry a geometry step seeds its next point with.
SavedState SnapshotState(const BasisSet::Complex_BS& bs, const WaveFunction::cWaveFunction& wf);
//! \c LowdinCarry over the solid's blocks.
void LowdinCarry(SavedState& carry, const BasisSet::Complex_BS& bs, SpinGroup g);

//! \brief LOAD a SOLID state: \c ReadState, refused (Kind) when the file holds another kind's state.
Outcome<SavedState, RestartRefusal> ReadSolidState(const std::string& path);

//! \brief Rebuild the saved density over THIS run's basis \a bs under \a g: one \c IrrepCD per block (the BZ weight re-applied), the
//! spin pairs projected (a magnetic group), composed exactly as the wave function composes its own.  The fingerprint must already
//! have been compared (the block sizes and k are ASSERTED against it here, not re-judged).  FAILS only on the Realness case.
Outcome<std::unique_ptr<ChargeDensity::cDM_CD>, RestartRefusal>
    RestoreDensity(const SavedState&, const BasisSet::Complex_BS& bs, SpinGroup g);

//  The file layout is qchem.SCFState's ("qchem-state 2.0", kind "solid"); each block carries its k[3] as a block extra.  A solid
//  "qchem-solid-state 1.x" file is read through the same reader.

} //export namespace qchem
