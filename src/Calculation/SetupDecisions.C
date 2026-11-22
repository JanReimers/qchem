// File: Calculation/SetupDecisions.C  The STRUCTURE-NEUTRAL decisions of a run setup, in ONE place every setup calls
// (doc/OOD-SOLID-Cleanup.md V-FACADE-SETUP).
//
// MoleculeSetup and SolidSetup each re-derived the same decisions in their own words, so a fix to one (the multiplicity-1 odd-count
// refusal, f59c7f65) had to be checked by hand against the other.  What is decided here reads PROPERTIES of the run -- an electron
// count, a multiplicity, whether the structure repeats, the run's pseudopotential set, an accelerator's name -- never the structure's
// type.  What stays in each setup is what is truly structure-specific: the k mesh, the cell and its strain box, the imposed group and
// its magnetic decoration, the XC-mesh sizing (an entry rung on a molecule, the sharpness gather on a cell), the basis recipe.
//
// Each decision is a policy of the policy layer: a refusal files its row (policy::RefuseAs) at the caller's report cursor, as before.
// The bodies live in Imp/SetupDecisions.C.
module;
#include <optional>
#include <string>
#include <vector>
#include <nlohmann/json.hpp>
export module qchem.SetupDecisions;
export import qchem.SeedPolicy;              // ResolveSeed, SeedKeyOf, SeedLibraryCovers: the seed's Auto rule (re-exported: one door)
import qchem.ChargeDensity.Seed;             // SeedStrategy
import qchem.ChargeDensity.AtomicDensity;    // SeedKey
import qchem.Structure;                      // Structure::MakeSymmetry (does the structure repeat?  HasLattice)
import qchem.Symmetry.Group;                 // HasLattice, kAdmissionTolerance
import qchem.Mesh;                           // qcMesh::MeshParams
import qchem.Pseudopotential.Set;            // PseudopotentialSet: the run's pseudo-atoms key the seed library
export import qchem.SCFAccelerator.Factory;  // SCFAccelerators::Type, AcceleratorParams (the accelerator block resolves to it)
export import qchem.SCFRecipe;               // SCFRecipe (AcceleratorOptions::recipe)
export import qchem.RunOptions;              // AcceleratorOptions (the block ResolveAcceleratorParams reads)
import qchem.Hamiltonian.Factory;            // AdmitVxc, VxcFit
export import qchem.Symmetry.Spin;           // SpinGroup: ResolveSpinGroup's answer

export namespace qchem
{

//=== SPIN BOOKKEEPING ===================================================================================================================
//! The electron count of each spin channel.
struct SpinChannels {int nUp = 0, nDown = 0;};

//! \brief THE MULTIPLICITY -> (nUp, nDown) RULE, every structure: 2S = nUp - nDown, nUp + nDown = \a nel.
//!   - \a multiplicity <= 0: MINIMAL spin -- 2S = nel mod 2 (a closed shell for an even count, a doublet for an odd one).
//!   - 1: the SINGLET, nUp = nDown, which an odd count cannot hold -- REFUSED there (until 2026-11-01 the solid mapped 1 to nel mod 2
//!     like 0, so the check could never fire: bcc Li, 1 e/cell, ran a fully spin-up doublet 53 mHa off and reported CONVERGED).
//!   - > 1: 2S = multiplicity - 1 as stated; REFUSED when 2S > nel or nel - 2S is odd (fail loud, never a silent miscount).
//! \a spinsShareFermi -- an option only where the structure has it (a solid's metal: ONE Fermi level sets the channels, so the counts
//! are only the split the occupations start from); \c nullopt where it is not one -- when true keeps the nel mod 2 split for a singlet
//! on an odd count instead of refusing it, and when false the refusal names it as the way out.  \a who prefixes the refusal (the options
//! struct the caller states).  Whether the run is then spin-POLARIZED is \c ResolveSpinGroup, below.
SpinChannels ResolveSpinChannels(int multiplicity, int nel, std::optional<bool> spinsShareFermi, const std::string& who);

//! \brief THE SPIN-GROUP RULE, every structure (user ruling R6, 2026-11-17, V-FACADE-PARITY): TWO ORTHOGONAL KNOBS.  The stated spin
//! group says whether the channels may differ (\c Polarized: two independent channels, so a local moment m(r) may be nonzero with zero
//! total moment -- an AFM crystal, a broken-symmetry dimer); the multiplicity says only S.  ONE rule ties them: a multiplicity > 1 is an
//! open shell, which an unpolarized run cannot hold, so it FORCES \c Polarized.  A broken-symmetry or AFM singlet, molecule or solid,
//! states \c Polarized with multiplicity 1.  (Until R6 the solid read multiplicity >= 1 as polarized and had no spin knob; a schema-1
//! deck's \c solid.multiplicity keeps that meaning through the deck's translator.)
SpinGroup ResolveSpinGroup(SpinGroup stated, int multiplicity);

//=== THE ACCELERATOR BLOCK ============================================================================================================
//! The accelerator BLOCK (\c AcceleratorOptions) lives in \c qchem.RunOptions; its resolution is here.

//! \brief \a p (the kind's defaults) with every knob \a acc STATES laid over it: the one place an accelerator block becomes concrete parameters.
SCFAccelerators::AcceleratorParams ResolveAcceleratorParams(const AcceleratorOptions& acc, SCFAccelerators::AcceleratorParams p);

//! \brief The MOLECULE's and the ATOM's defaults (the chain they always ran): Nproj 4, the absolute SV floor off, the relative prune
//! 1e-4, standalone GDM engaging at the same gate \a eMax as DIIS, trust 0.1, a Ladder GDM rung that always steps (gate 1e10).
SCFAccelerators::AcceleratorParams MolecularAcceleratorDefaults(double eMax);

//! The json door's keys for \a p (the atom's accelerator-sweep escape hatch merges over them; nothing else uses the json door).
nlohmann::json ToAcceleratorJson(const SCFAccelerators::AcceleratorParams& p);

//=== THE SEED ===========================================================================================================================
//! The run's resolved seed and the key of the seed library it reads.
struct RunSeed
{
    ChargeDensity::SeedStrategy seed = ChargeDensity::SeedStrategy::Default;
    ChargeDensity::SeedKey      key;
};
//! \brief THE SEED DECISION, every structure: the library key is the run's own pseudo-atoms' (\c SeedKeyOf; the unkeyed all-electron
//! library when \a pp is null), carrying the stated charge order and magnetic ordering, and the seed is \a stated or the Auto rule
//! (\c ResolveSeed) over the run's properties: does its iteration-0 Fock read a density (\a densityFock), does the structure repeat, are
//! formal charges stated, does the library hold every species.  Files the AutoSeed row when the rule runs.
RunSeed ResolveRunSeed(ChargeDensity::SeedStrategy stated, const Structure& st, const Pseudopotential::PseudopotentialSet* pp,
                              bool densityFock, const std::vector<int>& siteCharges = {}, const std::vector<int>& siteSpins = {});

//=== THE ACCELERATOR ====================================================================================================================
//! \brief The accelerator a run STATES BY NAME: DIIS | GDM | Ladder.  Any other name is REFUSED (row AcceleratorType), naming the three
//! -- never a silent DIIS.  \a who prefixes the refusal.
SCFAccelerators::Type AcceleratorTypeByName(const std::string& name, const std::string& who);
//! \brief THE MOLECULE'S AND THE ATOM'S ACCELERATOR, from its ONE home (\c RunOptions::accelerator, V-FACADE-PARITY P2): \c Auto is
//! their policy, \c DIIS; \c Null is REFUSED (policy AcceleratorType): with no Fock extrapolation a molecular SCF steps on the bare
//! diagonalization, and it has no density mixer to pair \c Null with, as the solid does.
SCFAccelerators::Type MolecularAcceleratorType(SCFAccelerators::Type asked, const std::string& who);
//! \brief DIIS's ENGAGE GATE (EMax, on [F,D]): \a stated when > 0, else the run's default.  A model whose Fock reads a density (DFT, a
//! pseudopotential run; \a densityFock) engages DIIS from iteration 0 (EMax 100): the Z-scaled gate keeps it off and the SCF limit-cycles
//! to a non-converged snapshot (the "M_Sym layout UB" mechanism, M_DFT).  Otherwise the proven Z-scaled gate \f$Z^2\cdot 0.08/32\f$
//! (0.08, was 0.1: the [F,D] it compares against changed UNIT on 2026-10-29 to the worst irrep, and on Kr the two differ by ~0.8 at
//! the engage step -- unrescaled, DIIS engaged one step early and the Kr libxc run went chaotic).
double DIISEngageGate(double stated, int Z, bool densityFock);

//=== THE Vxc ADMISSION (pin 30) =========================================================================================================
//! \brief Ask the Hamiltonian's admission rule (\c Hamiltonian::AdmitVxc) of the stated Vxc fit and XC mesh on \a st, BEFORE the basis
//! and the Hamiltonian are built: REFUSED (row AdmitVxc) with the rule's own reason.  \a where names the structure in the refusal.
void RequireVxcAdmitted(const Structure& st, Hamiltonian::VxcFit fit, const qcMesh::MeshParams& mesh, const std::string& where);

} //namespace qchem
