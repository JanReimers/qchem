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
module;
#include <optional>
#include <string>
#include <vector>
#include <nlohmann/json.hpp>
export module qchem.SetupDecisions;
export import qchem.SeedPolicy;              // ResolveSeed, SeedKeyOf, SeedLibraryCovers: the seed's Auto rule (re-exported: one door)
import qchem.ChargeDensity.Seed;             // SeedStrategy
import qchem.ChargeDensity.AtomicDensity;    // SeedKey
import qchem.Structure;                      // Structure::Admits (does the structure repeat?)
import qchem.Mesh;                           // qcMesh::ImposedSymmetry, MeshParams
import qchem.Pseudopotential.Set;            // PseudopotentialSet: the run's pseudo-atoms key the seed library
export import qchem.SCFAccelerator.Factory;  // SCFAccelerators::Type, AcceleratorParams (the accelerator block resolves to it)
export import qchem.SCFRecipe;               // SCFRecipe (AcceleratorOptions::recipe)
import qchem.Hamiltonian.Factory;            // AdmitVxc, VxcFit
import qchem.Policy;                         // policy::Refuse / RefuseAs
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
inline SpinChannels ResolveSpinChannels(int multiplicity, int nel, std::optional<bool> spinsShareFermi, const std::string& who)
{
    const bool singlet = multiplicity==1 && !spinsShareFermi.value_or(false);
    if (singlet && nel%2!=0)
        policy::RefuseAs("SpinMultiplicity", who+": multiplicity 1 is the SINGLET (nUp = nDown), impossible with an odd Nelec "
                                            +std::to_string(nel)+" -- state multiplicity 2 (a doublet) or 0 (minimal spin)"
                                            +(spinsShareFermi ? ", or spinsShareFermi true so one Fermi level sets the channels (a metal)" : ""));
    const int twoS = multiplicity>1 ? multiplicity-1 : (singlet ? 0 : nel%2);
    if (twoS>nel || (nel-twoS)%2!=0)
        policy::RefuseAs("SpinMultiplicity", who+": multiplicity "+std::to_string(multiplicity)+" (2S="+std::to_string(twoS)
                                            +") parity disagrees with Nelec "+std::to_string(nel)+" -- need 0 <= 2S <= Nelec and Nelec-2S even");
    return {(nel+twoS)/2, (nel-twoS)/2};
}

//! \brief THE SPIN-GROUP RULE, every structure (user ruling R6, 2026-11-17, V-FACADE-PARITY): TWO ORTHOGONAL KNOBS.  The stated spin
//! group says whether the channels may differ (\c Polarized: two independent channels, so a local moment m(r) may be nonzero with zero
//! total moment -- an AFM crystal, a broken-symmetry dimer); the multiplicity says only S.  ONE rule ties them: a multiplicity > 1 is an
//! open shell, which an unpolarized run cannot hold, so it FORCES \c Polarized.  A broken-symmetry or AFM singlet, molecule or solid,
//! states \c Polarized with multiplicity 1.  (Until R6 the solid read multiplicity >= 1 as polarized and had no spin knob; a schema-1
//! deck's \c solid.multiplicity keeps that meaning through the deck's translator.)
inline SpinGroup ResolveSpinGroup(SpinGroup stated, int multiplicity)
{
    return multiplicity>1 ? SpinGroup::Polarized : stated;
}

//=== THE ACCELERATOR BLOCK ============================================================================================================
//! \brief THE ACCELERATOR BLOCK of the options tree, every structure (V-FACADE-PARITY P0b-1, 2026-11-17; was the molecule's flat
//! \c AcceleratorOptions and the solid's \c SolidAcceleratorOptions).  Every knob is OPTIONAL: unstated takes the KIND's default, which
//! the kind's setup fills (\c ResolveAcceleratorParams) -- the molecule's and the atom's Z-scaled DIIS gate, Nproj 4 and a Ladder GDM rung
//! that always steps; the solid's GPW production recipe (\c SCFAccelerators::AcceleratorParams' member defaults).  Those differences
//! are per-kind DEFAULTS, logged for harmonizing (OpenWork "Per-kind defaults").  The knobs mirror \c DIISParams / \c GDMParams and
//! the Ladder's triggers one to one, so a knob added there has one spelling here.
struct AcceleratorOptions
{
    //! The accelerator BY NAME (molecule and atom; the solid's is \c SolidCalcOptions::accelerator until P0b-2): "DIIS" (default) |
    //! "GDM" (direct minimisation, robust for hard/open-shell SCF) | "Ladder" (DIIS -> GDM).  Any other name is REFUSED at Admit.
    std::string type = "DIIS";
    //! The SCF RECIPE by material regime (qchem.SCFRecipe; doc/SCFRobustnessPlan.md R4.1 as ruled for B5, ruling (c)): a molecule takes
    //! \c Auto or \c Insulator only; MagneticInsulator and Metal are REFUSED there (their Kerker+Pulay mixes in G space).
    SCFRecipe recipe = SCFRecipe::Auto;
    struct DIIS
    {
        std::optional<size_t> nProj;          //!< history depth
        //! DIIS starts extrapolating once [F,D] < this.  A molecule's / atom's unstated gate is Z-scaled (\c DIISEngageGate).
        std::optional<double> engageBelowFD;
        //! The ABSOLUTE floor on svMin(B) (OFF for a molecule since 2026-11-11: B scales as |[F,D]|^2, so it purged the history once
        //! [F,D] < ~7e-5) and the RELATIVE prune svMin(B) < svTolRel * max_i B_ii, which guards a genuinely dependent history (doc/Records/
        //! OpenWork_History5.md, "The molecular DIIS switches itself off").
        std::optional<double> svTol, svTolRel;
    };
    struct GDM
    {
        //! Engage the geodesic step once [F,D] < this.  A molecule's / atom's UNSTATED standalone-GDM gate is the DIIS gate (stated or
        //! Z-scaled), as the flat \c eMax always was; the Ladder's GDM rung keeps its own (1e10 there: it always steps).
        std::optional<double> engageBelowFD;
        std::optional<double> trust;          //!< the largest geodesic rotation (radians) per step
        std::optional<double> trustBackoff, trustMin, pcFloor, aufbauTol;
        std::optional<int>    aufbauMaxRepicks;
    };
    DIIS diis;
    GDM  gdm;      //!< standalone GDM
    //! LADDER (DIIS -> GDM): its GDM rung and the hand-off triggers (\c SCFAcceleratorLadder): energy still moving by more than
    //! \c eThresh for \c stall non-improving steps (STALL, \c floor its low backstop), or the signal below \c switchAt (TAIL; 0 = none;
    //! 1e-6 by default, the CLAUDE.md rule).  \c signal: |ΔE/E| (\c EnergyChange) or [F,D] (\c Error).
    struct Ladder
    {
        GDM gdm;
        std::optional<double> eThresh;
        std::optional<int>    stall;
        std::optional<double> floor, switchAt;
        std::optional<SCFAccelerators::ScheduleSignal> signal;
        bool startOnGDM = false;   //!< the seed is a RESTORED state: start on the GDM rung (a restart's own decision, never a user's)
    } ladder;
};


//! \brief \a kindDefaults with every knob \a acc STATES laid over it: the one place an accelerator block becomes concrete parameters.
//! \brief The MOLECULE's and the ATOM's defaults (the chain they always ran): Nproj 4, the absolute SV floor off, the relative prune
//! 1e-4, standalone GDM engaging at the same gate \a eMax as DIIS, trust 0.1, a Ladder GDM rung that always steps (gate 1e10).
//! The json door's keys for \a p (the atom's accelerator-sweep escape hatch merges over them; nothing else uses the json door).

inline SCFAccelerators::AcceleratorParams ResolveAcceleratorParams(const AcceleratorOptions& acc, SCFAccelerators::AcceleratorParams p)
{
    auto gdm=[](const AcceleratorOptions::GDM& g, SCFAccelerators::GDMParams& o)
    {
        if (g.engageBelowFD) o.EngageBelowFD=*g.engageBelowFD;
        if (g.trust)         o.Trust=*g.trust;
        if (g.trustBackoff)  o.TrustBackoff=*g.trustBackoff;
        if (g.trustMin)      o.TrustMin=*g.trustMin;
        if (g.pcFloor)       o.PCFloor=*g.pcFloor;
        if (g.aufbauTol)     o.AufbauTol=*g.aufbauTol;
        if (g.aufbauMaxRepicks) o.AufbauMaxRepicks=*g.aufbauMaxRepicks;
    };
    if (acc.diis.nProj)         p.diis.Nproj=*acc.diis.nProj;
    if (acc.diis.engageBelowFD) p.diis.EngageBelowFD=*acc.diis.engageBelowFD;
    if (acc.diis.svTol)         p.diis.SVTol=*acc.diis.svTol;
    if (acc.diis.svTolRel)      p.diis.SVTolRel=*acc.diis.svTolRel;
    gdm(acc.gdm, p.gdm);
    gdm(acc.ladder.gdm, p.ladderGdm);
    if (acc.ladder.eThresh)  p.ladderEThresh=*acc.ladder.eThresh;
    if (acc.ladder.stall)    p.ladderStall=*acc.ladder.stall;
    if (acc.ladder.floor)    p.ladderFloor=*acc.ladder.floor;
    if (acc.ladder.switchAt) p.ladderSwitchAt=*acc.ladder.switchAt;
    if (acc.ladder.signal)   p.ladderSignal=*acc.ladder.signal;
    p.ladderStartOnGDM = p.ladderStartOnGDM || acc.ladder.startOnGDM;
    return p;
}

inline SCFAccelerators::AcceleratorParams MolecularAcceleratorDefaults(double eMax)
{
    SCFAccelerators::AcceleratorParams p;
    p.diis      = {.Nproj=4, .EngageBelowFD=eMax, .SVTol=0.0, .SVTolRel=1e-4};
    p.gdm       = {.EngageBelowFD=eMax, .Trust=0.1};
    p.ladderGdm = {.EngageBelowFD=1e10, .Trust=0.1};    // the rung always steps once it is the active one
    p.ladderEThresh=1e-8; p.ladderStall=5; p.ladderFloor=1e-8; p.ladderSwitchAt=1e-6;
    p.ladderSignal=SCFAccelerators::ScheduleSignal::EnergyChange;
    return p;
}

inline nlohmann::json ToAcceleratorJson(const SCFAccelerators::AcceleratorParams& p)
{
    return {{"NProj",p.diis.Nproj},{"EngageBelowFD",p.diis.EngageBelowFD},{"SVTol",p.diis.SVTol},{"SVTolRel",p.diis.SVTolRel},{"Trust",p.ladderGdm.Trust},
            {"ethresh",p.ladderEThresh},{"stall",p.ladderStall},{"floor",p.ladderFloor},{"switchat",p.ladderSwitchAt},
            {"signal",p.ladderSignal==SCFAccelerators::ScheduleSignal::EnergyChange ? "EnergyChange" : "Error"}};
}

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
inline RunSeed ResolveRunSeed(ChargeDensity::SeedStrategy stated, const Structure& st, const Pseudopotential::PseudopotentialSet* pp,
                              bool densityFock, const std::vector<int>& siteCharges = {}, const std::vector<int>& siteSpins = {})
{
    RunSeed r;
    if (pp) r.key=SeedKeyOf(*pp, st);
    r.key.siteCharges=siteCharges;
    r.key.siteSpins  =siteSpins;
    r.seed=ResolveSeed(stated, st, {.densityFock = densityFock, .periodic = st.Admits(qcMesh::ImposedSymmetry::Translations),
                                    .chargesStated = !siteCharges.empty(),
                                    .covered = [key=r.key](int Z){ return SeedLibraryCovers(key, Z); }});
    return r;
}

//=== THE ACCELERATOR ====================================================================================================================
//! \brief The accelerator a run STATES BY NAME: DIIS | GDM | Ladder.  Any other name is REFUSED (row AcceleratorType), naming the three
//! -- never a silent DIIS.  \a who prefixes the refusal.
inline SCFAccelerators::Type AcceleratorTypeByName(const std::string& name, const std::string& who)
{
    using T=SCFAccelerators::Type;
    if (name=="DIIS")   return T::DIIS;
    if (name=="GDM")    return T::GDM;
    if (name=="Ladder") return T::Ladder;
    policy::RefuseAs("AcceleratorType", who+": unknown accelerator type \""+name+"\" (expected DIIS | GDM | Ladder)");
}
//! \brief DIIS's ENGAGE GATE (EMax, on [F,D]): \a stated when > 0, else the run's default.  A model whose Fock reads a density (DFT, a
//! pseudopotential run; \a densityFock) engages DIIS from iteration 0 (EMax 100): the Z-scaled gate keeps it off and the SCF limit-cycles
//! to a non-converged snapshot (the "M_Sym layout UB" mechanism, M_DFT).  Otherwise the proven Z-scaled gate \f$Z^2\cdot 0.08/32\f$
//! (0.08, was 0.1: the [F,D] it compares against changed UNIT on 2026-10-29 to the worst irrep, and on Kr the two differ by ~0.8 at
//! the engage step -- unrescaled, DIIS engaged one step early and the Kr libxc run went chaotic).
inline double DIISEngageGate(double stated, int Z, bool densityFock)
{
    return stated>0.0 ? stated : (densityFock ? 100.0 : Z*Z*0.08/32);
}

//=== THE Vxc ADMISSION (pin 30) =========================================================================================================
//! \brief Ask the Hamiltonian's admission rule (\c Hamiltonian::AdmitVxc) of the stated Vxc fit and XC mesh on \a st, BEFORE the basis
//! and the Hamiltonian are built: REFUSED (row AdmitVxc) with the rule's own reason.  \a where names the structure in the refusal.
inline void RequireVxcAdmitted(const Structure& st, Hamiltonian::VxcFit fit, const qcMesh::MeshParams& mesh, const std::string& where)
{
    if (const auto a=Hamiltonian::AdmitVxc(st, fit, mesh); !a)
        policy::Refuse("AdmitVxc", "the stated Vxc fit and XC mesh on this "+where, a.Error());
}

} //namespace qchem
