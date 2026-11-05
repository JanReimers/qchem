// File: Calculation/SolidState.C  A periodic SCF state ON DISK: what is saved, the fingerprint that decides
// whether a later run may start from it, and the rebuild of the density it carries (doc/OpenWork.md §2
// "SCF checkpoint/restart", CK-1).
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
#include <cstdint>
#include <limits>
#include <map>
#include <memory>
#include <string>
#include <vector>
export module qchem.SolidState;
export import qchem.Outcome;
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

//! \brief WHY a saved state may not seed this run.  A value, not a throw: a campaign that finds no usable
//! state falls back to a fresh seed, which is a decision only the caller can make.
struct RestartRefusal
{
    enum class Why
    {
        Unreadable,   //!< absent, not HDF5, or not a qchem solid state
        Format,       //!< a qchem state of a format version this build does not read
        Mismatch,     //!< a REFUSE-class fingerprint entry differs (details name every one)
        Realness      //!< a COMPLEX saved block restarting into a REAL block whose D is not real
    };
    Why         why = Why::Unreadable;
    std::string details;
};

//! \brief What identifies a run to a saved state, in the two classes the file header explains.
//!
//! Plain named arrays rather than a struct of typed fields: each entry is written as its own dataset, so a
//! Python reader sees them by name, and a new entry is one line in \c MakeStateFingerprint -- the comparison,
//! the writer and the reader are generic over the maps.
struct StateFingerprint
{
    std::map<std::string, std::vector<double>> refuseNum, warmNum;
    std::map<std::string, std::string>         refuseStr, warmStr;
};

//! The inputs of a fingerprint that only the FACADE knows (its options); the rest is asked of the built basis.
struct RunIdentity
{
    std::vector<std::pair<std::string,int>>      species;      //!< (element, valence)
    int                                          Nelec = 0;
    int                                          multiplicity = 0;
    std::string                                  functional = "LDA";   //!< recorded, WARM (an LDA->PBE restart is a legitimate continuation)
    std::string                                  ppSet;        //!< the pseudopotential SET's name (the family; q per Z is in \c species)
    std::vector<Hamiltonian::HubbardManifold>    hubbard;      //!< the U the Hamiltonian CURRENTLY carries
    double                                       densityEcut = 0.0, cutoffFactor = 0.0;
    std::string                                  xcMesh;       //!< the RESOLVED XC quadrature, named
    int                                          ortho = 0;
    double                                       orthoTol = 0.0;
    double                                       kT = 0.0;     //!< the recipe's temperature (warm)
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

//! The fingerprint of a run: the cell matrix \a A (columns = lattice vectors), the sites of \a st, the built
//! Bloch basis \a bs (its blocks' k, weight and size; the orbital shells off the AoShellSource face) under spin
//! group \a g, plus the facade's \a id.
StateFingerprint MakeStateFingerprint(const rmat3d_t& A, const Structure& st, const BasisSet::Complex_BS& bs,
                                      SpinGroup g, const RunIdentity& id);

//! \brief The outcome of comparing a saved fingerprint against this run's: REFUSED on any refuse-class
//! difference, else the warm-class differences (EMPTY = an exact resume), one readable line each.
Outcome<std::vector<std::string>, RestartRefusal> CompareFingerprints(const StateFingerprint& saved,
                                                                     const StateFingerprint& now);

//! What the run reports about itself beside the state (echoed into the file, read back on restart).
struct StateSummary
{
    std::string label;
    bool        converged  = false;
    size_t      iterations = 0;
    double      energy     = 0.0;     //!< total energy, hartree (the LAST iterate's when not converged)
    double      charge     = 0.0;
    double      commutator = 0.0;     //!< the final iteration's [F,D] (NaN = unmeasured: a Null accelerator, i.e. a density-MIXED run)
    std::string recipe;               //!< the last stage's SCF banner facts, human-readable
    //! The final iteration's relative density change -- the residual the state was SAVED at (format 1.1, 2026-10-30).  NaN in a
    //! format-1.0 state, which did not record it: an OPTIONAL payload field, the first one (CleanCode D-STATE-SCHEMA).
    double      densityChange = std::numeric_limits<double>::quiet_NaN();
};

//! \brief SAVE: the fingerprint, the summary, and every (k,σ) block's D, C, ε, f from \a wf over \a bs.
//! Written to \a path + ".tmp", flushed, then RENAMED into place -- a crash mid-write never destroys the
//! previous good state.  THROWS on an I/O failure.
void WriteSolidState(const std::string& path, const StateFingerprint&, const StateSummary&,
                     const BasisSet::Complex_BS& bs, const WaveFunction::cWaveFunction& wf);

//! One saved block, as read back.  D is always held complex here (a real block's D promotes exactly).
struct SavedBlock
{
    std::int64_t          ms = 0;         //!< the Spin, as its enum value
    rvec3_t               k;              //!< fractional
    double                weight = 0.0;   //!< the block's BZ weight
    bool                  real = false;   //!< the SAVING run's working type for this block
    mat_t<dcmplx>         D;              //!< physical (unweighted) AO density matrix
    //! The SAVING basis's overlap -- IN MEMORY ONLY (\c SnapshotState fills it, a file never carries it): what \c LowdinCarry
    //! needs to take D from one geometry's basis to the next.  Empty after \c ReadSolidState.
    mat_t<dcmplx>         S;
};

//! A state read from disk.
struct SavedState
{
    std::string             path;
    StateFingerprint        fingerprint;
    StateSummary            summary;
    std::vector<SavedBlock> blocks;       //!< in the canonical order (spin irreps outer, basis blocks inner)
};

//! \brief The IN-MEMORY half of \c WriteSolidState: every (k,σ) block's D from \a wf over \a bs, in the canonical block order,
//! exactly as \c RestoreDensity reads it -- the CARRY a geometry step seeds its next point with (a restart from memory, whose
//! basis is the same recipe at moved atoms).  The fingerprint and summary are left empty: the caller built both points from one
//! setup, which is the comparison a file's fingerprint stands in for.  Each block also keeps its overlap S (\c LowdinCarry).
SavedState SnapshotState(const BasisSet::Complex_BS& bs, const WaveFunction::cWaveFunction& wf);

//! \brief THE LÖWDIN CARRY (doc/CalcNetworkReview.md §5 step 8): take a snapshot's D from the basis it was made on to \a bs -- the SAME
//! recipe at moved atoms -- through the Löwdin-orthonormal frame, \f$D' = A D A^\dagger,\ A = S_{new}^{-1/2} S_{old}^{1/2}\f$, block for
//! block.  The orbitals' Löwdin coefficients \f$S^{1/2}C\f$ are held fixed, so D' stays idempotent in the NEW metric
//! (\f$D'S'D' = D'\f$, \f$\mathrm{Tr}\,D'S' = N\f$), which the verbatim AO carry is not once the centres move.  THROWS if a block
//! has no saved overlap (a state read from a file) or its size differs from \a bs's.  The block walk is \c SnapshotState's.
void LowdinCarry(SavedState& carry, const BasisSet::Complex_BS& bs, SpinGroup g);

//! \brief LOAD.  FAILS on an absent/foreign file or an unreadable format; THROWS on a file that claims to be
//! a qchem state and is internally inconsistent (a broken invariant: somebody edited it).
Outcome<SavedState, RestartRefusal> ReadSolidState(const std::string& path);

//! \brief Rebuild the saved density over THIS run's basis \a bs under \a g: one \c IrrepCD per block (the
//! BZ weight re-applied), composed exactly as the wave function composes its own.  The fingerprint must
//! already have been compared (the block sizes and k are ASSERTED against it here, not re-judged).  FAILS
//! only on the Realness case.
Outcome<std::unique_ptr<ChargeDensity::cDM_CD>, RestartRefusal>
    RestoreDensity(const SavedState&, const BasisSet::Complex_BS& bs, SpinGroup g);

//  THE FILE LAYOUT (format "qchem-solid-state 1.1"; the reader takes any 1.x -- a minor version only ADDS optional fields):
//    /                     attrs format, created (UTC), label, converged, iterations, energy, charge,
//                                commutator, recipe, densityChange (1.1; absent in 1.0)
//    /fingerprint/refuse   one dataset per numeric entry; string entries as attributes
//    /fingerprint/warm     the same
//    /blocks/<b>           attrs ms, weight, real, n, nmo; datasets k[3], D[n,n], C[n,nmo], eps[nmo], f[nmo]
//                          (C order; complex as the {r,i} compound; f physical occupations)

} //export namespace qchem
