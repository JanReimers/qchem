// File: Calculation/SeedPolicy.C  The SCF seed's Auto rule: ONE function for every structure (user, 2026-11-12: "we should strive to make
// seeding and spherical/Cartesian options structure neutral"; rulings of the same day: "(1) Derived only if periodic", "Core guess + WARN,
// both").
//
// The rule, in order (a stated seed is honoured and never reaches it):
//   1. A model whose iteration-0 Fock reads no density (Hartree-Fock / one-electron, all-electron): the core-Hamiltonian guess.
//   2. A species the run's seed library lacks: the core guess, and a WARN -- the seed the rule would have chosen cannot be built.
//   3. Formal charges STATED: the SAD superposition at those charges (IonicSAD).
//   4. A PERIODIC structure: the SAD superposition at the charges the electronegativity heuristic derives (IonicSAD).  WHY PERIODICITY:
//      a lattice's Madelung field stabilises the ionic charge states (Na+ F- in rock salt), so an ionic start is the better guess there
//      (V2.2: Uniform has a stable wrong basin for electron-sparse cells, the neutral SAD leaves the whole charge transfer to the SCF).  An
//      isolated molecule has no such field: MEASURED 2026-11-12 (LDA/dzvp, M_DFT.*_IonicSADReachesTheSADState), the derived charges reach
//      the same state but cost iterations even on an ionic molecule -- H2O 8 -> 10, LiF 8 -> 10, OH- 9 -> 10.
//   5. Otherwise: the neutral SAD superposition (SAD).
// The structure's TYPE never enters: the rule reads properties (does it repeat, which species, are charges stated), and both setups call
// it.  A policy of the policy layer: it files one `policy` row (AutoSeed) at the caller's report cursor.
module;
#include <functional>
#include <iostream>
#include <string>
export module qchem.SeedPolicy;
import qchem.ChargeDensity.Seed;   // SeedStrategy
import qchem.Structure;            // Structure (its species)
import qchem.PeriodicTable;        // GetSymbol (the warning names the species)
import qchem.Policy;
import qchem.Pseudopotential.Set;  // PseudopotentialSet: the run's pseudo-atoms, which key the seed library

export namespace qchem
{

//! The seed library key for a run, the same for every structure: the set's functional and, for each element of \a st, the GTH variant q
//! its Hamiltonian uses -- so the seed density is the SAME pseudo-atom as the terms (C4 step 3; for a molecule, OpenWork §3 "A
//! pseudopotential run seeds from the VALENCE library for its q, on every structure", 2026-11-12).  An all-electron run has no set and
//! uses the unkeyed \c SeedKey{} (the all-electron library, \c ChargeDensity::SeedLibrary).
inline ChargeDensity::SeedKey SeedKeyOf(const Pseudopotential::PseudopotentialSet& pp, const Structure& st)
{
    ChargeDensity::SeedKey key(pp.Functional());
    for (const auto& [element, q] : pp.Species(st)) key.q[(int)thePeriodicTable().GetZ(element)] = q;
    return key;
}
//! Does the run's seed library (\c ChargeDensity::SeedLibrary) hold species \a Z?  The Auto rule's coverage, both setups.
inline bool SeedLibraryCovers(const ChargeDensity::SeedKey& key, int Z)
{
    return ChargeDensity::HasAtomicDensity(Z, key.functional, ChargeDensity::SeedLibrary(key,Z), -1, key.Q(Z));
}

//! What the Auto rule reads about a run.  Properties, never the structure's type.
struct SeedFacts
{
    bool densityFock;                     //!< the iteration-0 Fock reads a density (DFT, or any pseudopotential run)
    bool periodic;                        //!< the structure is a lattice (its Madelung field holds ionic charge states)
    bool chargesStated;                   //!< the run states each site's formal charge
    std::function<bool(int Z)> covered;   //!< the run's seed library holds species \a Z
};

//! The run's seed: \a stated when it is not \c Default, else the Auto rule above, which files its row (and prints a WARN to the console).
inline ChargeDensity::SeedStrategy ResolveSeed(ChargeDensity::SeedStrategy stated, const Structure& st, const SeedFacts& f)
{
    using ChargeDensity::SeedStrategy;
    if (stated!=SeedStrategy::Default) return stated;
    if (!f.densityFock)
    {
        policy::Report("AutoSeed", "seed: the core-Hamiltonian guess", policy::Verdict::Choose,
                       "a Hartree-Fock / one-electron run, and no seed stated");
        return SeedStrategy::CoreGuess;
    }
    std::string missing;
    st.ForEachSite([&](int Z, const rvec3_t&, bool)
    {
        const std::string sym=thePeriodicTable().GetSymbol(Z);
        if (!f.covered(Z) && missing.find(" "+sym+" ")==std::string::npos) missing+=" "+sym+" ";
    });
    if (!missing.empty())
    {
        const std::string decided="seed: the core-Hamiltonian guess (the seed library lacks"+missing+")";
        const std::string reason="no seed stated, and the SAD superposition the run would take needs every species in its seed library; "
                                 "generate the entry (CLIapps/valgen) or state a seed";
        policy::Report("AutoSeed", decided, policy::Verdict::Warn, reason);
        std::cout << "[policy] WARNING " << decided << " -- " << reason << std::endl;
        return SeedStrategy::CoreGuess;
    }
    if (f.chargesStated)
    {
        policy::Report("AutoSeed", "seed: IonicSAD, the superposition at the stated formal charges", policy::Verdict::Choose,
                       "formal charges stated, and no seed stated");
        return SeedStrategy::IonicSAD;
    }
    if (f.periodic)
    {
        policy::Report("AutoSeed", "seed: IonicSAD, the superposition at the derived formal charges", policy::Verdict::Choose,
                       "a periodic structure: its Madelung field holds ionic charge states, so the derived charges are the better start");
        return SeedStrategy::IonicSAD;
    }
    policy::Report("AutoSeed", "seed: SAD, the superposition of neutral-atom densities", policy::Verdict::Choose,
                   "an isolated structure, no charges stated: no Madelung field holds an ionic start (it costs iterations, measured)");
    return SeedStrategy::SAD;
}

} //namespace
