// File: Calculation/SCFState.C  An SCF state ON DISK, for every kind of calculation (atom, molecule, solid): what is saved, the
// fingerprint that decides whether a later run may start from it, the file framing, and the carry of a density across a moved basis.
//
// V-FACADE-STATE (doc/OOD-SOLID-Cleanup.md; user ruling via lane G 2026-11-02, executed by lane A 2026-11-09): the abstract
// \c Calculation saves and restores what every kind SHARES -- the wave-function blocks (D, C, ε, f per Irrep), the run summary, the
// basis identity, the framing -- and each facade adds its own fingerprint entries.  This module is that shared half, moved out of
// \c qchem.SolidState (which keeps only the periodic parts: the k-blocks, the cell, the GPW knobs, the Shubnikov spin-pair restore).
//
// WHAT IS SAVED, per block (spin irreps outer, basis blocks inner -- the order a composite wave function builds its children in):
//   D     the block's AO density matrix, PHYSICAL occupations (a solid's BZ weight is stored beside it, never in it)
//   C,ε,f ALL orbitals (occupied and virtual): a WaveFunction read from disk (χ₀ / ACBN0 / gaps on a stored state with no SCF)
//   S     the block's overlap (2.0): what the Löwdin carry needs to take D onto a basis at MOVED atoms (\c LowdinCarry)
// plus the FINGERPRINT and the run's own summary.  No mixer or accelerator history: a restart is a FRESH STAGE (the 2026-09-21
// stale-Pulay lesson).
//
// THE FINGERPRINT has two classes, and the split IS the design:
//   REFUSE  the saved D is a matrix over THESE functions for THIS many electrons: the atom kinds, the orbital basis, the blocks,
//           Nelec, the spin group.  Any difference and the state means nothing here -- a REFUSED Outcome, never a wrong start.
//   WARM    a restart ACROSS these is what the mechanism is FOR: the functional, the multiplicity, the meshes, the atom POSITIONS (a
//           geometry step: the density is carried through the Löwdin frame, so Tr(D S') = N on the moved basis).  No difference at
//           all = an EXACT RESUME, which starts AT the saved answer.
module;
#include <cstdint>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <string>
#include <vector>
export module qchem.SCFState;
export import qchem.Outcome;
import qchem.Types;
import qchem.Structure;
import qchem.Symmetry.Irrep;             // Irrep (a block's spin and symmetry)
import qchem.Symmetry.Spin;              // SpinGroup
import qchem.BasisSet.Orbital_1E_IBS;    // the per-block bases (real and complex)
import qchem.Orbitals;                   // Orbitals (a block's C, ε, f)

export namespace qchem
{

//! \brief WHY a saved state may not seed this run.  A value, not a throw: a campaign that finds no usable state falls back to a
//! fresh seed, which is a decision only the caller can make.
struct RestartRefusal
{
    enum class Why
    {
        Unreadable,   //!< absent, not HDF5, or not a qchem state
        Format,       //!< a qchem state of a format version this build does not read
        Kind,         //!< a state of another KIND of calculation (a molecule's state offered to an atom run)
        Mismatch,     //!< a REFUSE-class fingerprint entry differs (details name every one)
        Realness      //!< a COMPLEX saved block restarting into a REAL block whose D is not real
    };
    Why         why = Why::Unreadable;
    std::string details;
};

//! \brief What identifies a run to a saved state, in the two classes the file header explains.  Plain named arrays: each entry is
//! its own dataset, so a Python reader sees them by name, and a new entry is one line in a facade's fingerprint.
struct StateFingerprint
{
    std::map<std::string, std::vector<double>> refuseNum, warmNum;
    std::map<std::string, std::string>         refuseStr, warmStr;
};

//! What the run reports about itself beside the state (echoed into the file, read back on restart).
struct StateSummary
{
    std::string label;
    bool        converged  = false;
    size_t      iterations = 0;
    double      energy     = 0.0;     //!< total energy, hartree (the LAST iterate's when not converged)
    double      charge     = 0.0;
    double      commutator = std::numeric_limits<double>::quiet_NaN();   //!< the final [F,D] (NaN = unmeasured: a density-MIXED run)
    std::string recipe;               //!< the last stage's SCF banner facts, human-readable
    //! The final iteration's relative density change -- the residual the state was SAVED at.  NaN in a state that did not record it.
    double      densityChange = std::numeric_limits<double>::quiet_NaN();
};

//! \brief ONE BLOCK OF THE RUN, as the facade hands it to the shared half: its Irrep, its basis (exactly one of \c real / \c cplx),
//! the orbitals the wave function holds for it, and the per-block identity only its kind knows (\c extras: a solid's k).
struct StateBlock
{
    Irrep                                    irrep;
    const BasisSet::Orbital_1E_IBS<double>*  real     = nullptr;   //!< set on a REAL block, else null
    const BasisSet::Orbital_1E_IBS<dcmplx>*  cplx     = nullptr;   //!< set on a complex block, else null
    const Orbitals::Orbitals*                orbitals = nullptr;   //!< null on a walk that needs no orbitals (a fingerprint, a restore)
    double                                   weight   = 1.0;       //!< the block's BZ weight (1 for a molecule or an atom)
    std::map<std::string, std::vector<double>> extras;             //!< written as datasets of the block, read back into \c SavedBlock
    size_t n() const;                                              //!< the block's function count
    std::string Label() const;                                     //!< the spatial irrep's label (the spin is \c irrep.ms)
};

//! One saved block, as read back.  D is always held complex here (a real block's D promotes exactly).
struct SavedBlock
{
    std::int64_t          ms = 0;         //!< the Spin, as its enum value
    std::string           irrep;          //!< the spatial irrep's label (empty in a solid 1.x state, which did not record it)
    double                weight = 1.0;   //!< the block's BZ weight
    bool                  real = false;   //!< the SAVING run's working type for this block
    std::map<std::string, std::vector<double>> extras;   //!< a solid's "k" (fractional)
    mat_t<dcmplx>         D;              //!< physical (unweighted) AO density matrix
    //! The SAVING basis's overlap: in a 2.0 file, and in an in-memory snapshot.  EMPTY from a solid 1.x file -- \c LowdinCarry
    //! cannot move such a block.
    mat_t<dcmplx>         S;
};

//! A state read from disk (or snapshot in memory).
struct SavedState
{
    std::string             path;
    std::string             kind;         //!< "atom" | "molecule" | "solid"
    StateFingerprint        fingerprint;
    StateSummary            summary;
    std::vector<SavedBlock> blocks;       //!< in the canonical order (spin irreps outer, basis blocks inner)
};

//! \name THE BASE FINGERPRINT: what every kind records the same way (the keys a solid state has always carried)
//!@{
//! REFUSE structure.Z; WARM structure.positions.
void AddStructure(StateFingerprint&, const Structure&);
//! REFUSE electrons.Nelec, spinGroup; WARM functional, spin.multiplicity.
void AddElectrons(StateFingerprint&, int Nelec, SpinGroup, int multiplicity, const std::string& functional);
//! REFUSE blocks.ms, blocks.n, blocks.irrep (the spatial labels, joined); WARM blocks.real (a working-type change promotes D exactly).
void AddBlocks(StateFingerprint&, const std::vector<StateBlock>&);
//! The orbital SHELLS of an atom-centred Gaussian basis (\c AoShellSource), in the order given: REFUSE basis.shellL / Components /
//! Primitives / exponents / coefficients, WARM basis.shellCenters (they follow the atoms).  THROWS if \a blocks is empty or a
//! block does not report its shells.
void AddShells(StateFingerprint&, const std::vector<const BasisSet::Orbital_1E_IBS<double>*>& blocks);
void AddShells(StateFingerprint&, const BasisSet::Orbital_1E_IBS<dcmplx>& block);
//!@}

//! \brief A migration of a saved fingerprint written under an older schema of one KIND, given this run's (\a now), applied before the
//! comparison.  Additive: an old file keeps loading.  Each kind owns its own (the solid's: \c SolidState's \c MigrateSolidFingerprint).
using FingerprintMigration = std::function<StateFingerprint(StateFingerprint saved, const StateFingerprint& now)>;

//! \brief The outcome of comparing a saved fingerprint against this run's: REFUSED on any refuse-class difference, else the
//! warm-class differences (EMPTY = an exact resume), one readable line each, each with what it MEANS for the restart.  A moved atom
//! is one line (the max displacement), not one per coordinate.  \a migrate, when given, is applied to \a saved first.
Outcome<std::vector<std::string>, RestartRefusal> CompareFingerprints(const StateFingerprint& saved, const StateFingerprint& now,
                                                                     const FingerprintMigration& migrate = nullptr);
//! The largest atom displacement between the two fingerprints' positions (bohr), 0 when either has none.
double MaxDisplacement(const StateFingerprint& saved, const StateFingerprint& now);

//! \brief SAVE: the framing (format, kind, the summary), the fingerprint, and every block's D, C, ε, f, S (every \a blocks entry
//! needs its orbitals).  Written to \a path + ".tmp", flushed, then RENAMED into place -- a crash mid-write never destroys the
//! previous good state.  THROWS on an I/O failure.
void WriteState(const std::string& path, const std::string& kind, const StateFingerprint&, const StateSummary&,
                const std::vector<StateBlock>& blocks);

//! \brief LOAD.  Reads every "qchem-state 2.x", and a "qchem-solid-state 1.x" as kind "solid" (its block k into \c extras, no
//! irrep labels, no S).  FAILS on an absent/foreign file or an unreadable format; THROWS on a file that claims to be a qchem state
//! and is internally inconsistent (a broken invariant: somebody edited it).
Outcome<SavedState, RestartRefusal> ReadState(const std::string& path);

//! \brief The IN-MEMORY half of \c WriteState: every block's D and S, exactly as a restore reads them.  The fingerprint and summary
//! are left empty (the caller built both points from one setup, which is the comparison a file's fingerprint stands in for).
SavedState SnapshotState(const std::vector<StateBlock>& blocks, const std::string& kind);

//! \brief THE LÖWDIN CARRY: take each block's D from the basis it was saved on to \a blocks' -- the SAME recipe at moved atoms --
//! through the Löwdin-orthonormal frame, \f$D' = A D A^\dagger,\ A = S_{new}^{-1/2} S_{old}^{1/2}\f$.  The orbitals' Löwdin
//! coefficients \f$S^{1/2}C\f$ are held fixed, so \f$\mathrm{Tr}\,D'S' = \mathrm{Tr}\,DS = N\f$ for ANY occupations, and D' stays
//! idempotent in the new metric where D was in the old -- which the verbatim AO carry is not once the centres move.  THROWS if a
//! block has no saved overlap (a solid 1.x file) or its size differs.
void LowdinCarry(SavedState& carry, const std::vector<StateBlock>& blocks);

//! \brief WHAT A RESTART DID, for the run's clients (a GUI shows it; a gate checks it): where it came from, whether it resumed
//! exactly, every warm difference (one line each, with what it means), how far the atoms moved, and the electrons the restored
//! density holds on THIS run's basis, \f$\sum_b w_b\,\mathrm{Re\,Tr}(D_b S'_b)\f$ -- N when the carry is right.
struct RestartSummary
{
    std::string              path;
    bool                     exact = false;
    std::vector<std::string> warm;                  //!< empty = an exact resume
    double                   maxDisplacement = 0.0; //!< bohr
    bool                     lowdinCarried = false; //!< the density was carried through the Löwdin frame (moved atoms, S saved)
    double                   electrons = 0.0;       //!< Σ_b w_b Re Tr(D_b S'_b) of the restored density, on this run's basis
};

//! \f$\sum_b w_b\,\mathrm{Re\,Tr}(D_b S_b)\f$ of \a st's densities over \a blocks' overlaps (the block weights are \a blocks').
double Electrons(const SavedState& st, const std::vector<StateBlock>& blocks);

//! \brief THE RESTORE CHECKS every kind shares, on a state whose fingerprint has been compared: one saved block per run block, with
//! its size and spin (a disagreement here is a fingerprint missing an entry: a THROW), and the REALNESS of a block restarting real
//! from a complex save (a REFUSAL: |Im D| > 1e-8).  Hands back each block's D, hermitised and promoted to complex.
Outcome<std::vector<mat_t<dcmplx>>, RestartRefusal> CheckedBlockDensities(const SavedState&, const std::vector<StateBlock>& blocks);

//  THE FILE LAYOUT (format "qchem-state 2.0"; the reader takes any 2.x, and the solid's 1.x):
//    /                     attrs format, kind, created (UTC), label, converged, iterations, energy, charge, commutator, recipe,
//                                densityChange
//    /fingerprint/refuse   one dataset per numeric entry; string entries as attributes
//    /fingerprint/warm     the same
//    /blocks/<b>           attrs ms, irrep, weight, real, n, nmo; datasets D[n,n], C[n,nmo], eps[nmo], f[nmo], S[n,n], and one per
//                          block extra (a solid block: k[3])
//                          (C order; complex as the {r,i} compound; f physical occupations)
//  The solid's 1.x layout is this one without kind, irrep and S (and format "qchem-solid-state 1.x").

} //export namespace qchem
