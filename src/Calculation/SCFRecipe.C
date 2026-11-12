// File: Calculation/SCFRecipe.C  THE SCF RECIPE POLICY: which accelerator and which mixer knobs a stage runs, from the MATERIAL REGIME
// the caller names (or Auto) and the STAGE KIND the facade knows.  doc/SCFRobustnessPlan.md §3 (the regime matrix) and R4.1 as ruled
// for B5 (user, 2026-11-12).
//
// A recipe names the REGIME only (Insulator | MagneticInsulator | Metal, or Auto).  The stage kind (cold, carried, restored) is never a
// name a user picks: the facade always knows it, and every regime has its cold / carried / restored column (ruling (a)).  So a recipe
// resolves through a 2-D table, (regime) x (stage kind), and that table IS §3:
//
//   regime \ stage      COLD (seeded)                  CARRIED (anneal, U loop, +-alpha, geometry)   RESTORED (same Hamiltonian, row 4)
//   Insulator (row 1)   Null + Pulay-D (0.45, 8, 5)    Null + Pulay-D, fresh history (row 5)          kT 0: Ladder@GDM; kT > 0 or a partial
//   MagneticInsulator   Null + K+P (G0 1, 0.45, 8, 5)  Null + K+P, fresh history                        eigensolve: DIIS.  Either beside LINEAR
//     (row 2)                                                                                           D-mixing, its START <= 0.3
//   Metal (row 3)       Null + K+P (G0 1, PROVISIONAL) Null + K+P, fresh history
//
// AUTO, in B5, IS TODAY'S POLICY NAMED AND PRINTED, and nothing else: it resolves every stage to the Insulator column, which is exactly what
// SolidCalculation::MakeAccelerator + qchem::ResolveMixing did before (the golden table in tests/SCFRecipe.C holds the old rules verbatim).
// Auto picking Metal (after Q4's Al measurement) or MagneticInsulator (after the C1-DEG detector) is LATER work that moves anchors.
//
// SMEARING IS A MODEL CHOICE, so no recipe sets kT: MagneticInsulator and Metal REFUSE kT = 0 (§3 row 2', ruling (b)); Insulator takes any.
// STATED KNOBS WIN: a knob the caller set (not the SCFParams default) beside a named recipe is kept and printed as a DEVIATION, so every
// existing deck and test means exactly what it meant (decks stay backward compatible: user convention 2026-11-12).
module;
#include <algorithm>
#include <iostream>
#include <optional>
#include <sstream>
#include <string>
#include <vector>
export module qchem.SCFRecipe;
export import qchem.SCFParams;
export import qchem.SCFAccelerator.Factory;   // SCFAccelerators::Type
import qchem.Policy;
import qchem.Reporting;

export namespace qchem
{

//! The MATERIAL REGIME a stage's SCF recipe is chosen for (doc/SCFRobustnessPlan.md §3; deck key `scf.recipe`).
enum class SCFRecipe { Auto, Insulator, MagneticInsulator, Metal };
//! What the stage starts from -- the facade's knowledge, never the user's name (R4.1 ruling (a)).
enum class SCFStageKind { Cold, Carried, Restored };

inline const char* NameOf(SCFRecipe r)
{
    switch (r)
    {
    case SCFRecipe::Auto:              return "Auto";
    case SCFRecipe::Insulator:         return "Insulator";
    case SCFRecipe::MagneticInsulator: return "MagneticInsulator";
    case SCFRecipe::Metal:             return "Metal";
    }
    return "?";
}
inline const char* NameOf(SCFAccelerators::Type t)
{
    using T=SCFAccelerators::Type;
    switch (t)
    {
    case T::DIIS:   return "DIIS";
    case T::GDM:    return "GDM";
    case T::Ladder: return "Ladder";
    case T::Null:   return "Null";
    case T::Auto:   return "Auto";
    }
    return "?";
}
inline const char* NameOf(SCFStageKind k)
{
    switch (k)
    {
    case SCFStageKind::Cold:     return "COLD";
    case SCFStageKind::Carried:  return "CARRIED";
    case SCFStageKind::Restored: return "RESTORED";
    }
    return "?";
}

//! Rows 2 and 3's mixer: Kerker-preconditioned Pulay on rho-tilde -- exactly the banked free MnO AFM-II recipe's knobs
//! (`decks/MnO_AFM2_free.json`, CLAUDE.md; doc/Benchmark.md the CP2K-LOOP rows).  ⚠ Metal's G0 1.0 is PROVISIONAL: the K+P decks'
//! value, UNMEASURED on a large metallic cell (§3 row 3) until Q4's Al measurement (Pulay-D vs K+P).
inline constexpr MixingRecipe theKerkerPulayRecipe{.relax=0.45, .G0=1.0, .depth=8, .start=5};
//! Row 4's linear step on a restored density: a converged density wants a restart's step, not the ground state's (relax 1.0 diverged a
//! restart, Si Γ 2026-09-29; CLAUDE.md "a linear fallback relaxes at 0.2-0.3").  ⚠ Applied to the START only, as before B5; OpenWork's
//! NEAR-TERM row "unify the restored stage's linear cap on the CEILING" lands that as its own measured step.
inline constexpr double theRestoredLinearStep=0.3;

//! The regime the policy RUNS: Auto resolves to Insulator in B5 (see the file comment).
inline SCFRecipe Resolved(SCFRecipe r) {return r==SCFRecipe::Auto ? SCFRecipe::Insulator : r;}

//! \brief REFUSE a recipe the stage cannot run as asked (a policy::Refuse: its row, then std::invalid_argument).
//! \a molecule: MagneticInsulator / Metal have no molecular column (no G-space Kerker on a molecule; ruling (c)).  At kT = 0 they refuse:
//! a near-degenerate d frontier has no integer aufbau solution (§3 row 2': J5, F4), and a metal's sloshing needs the smeared fill.
inline void AdmitRecipe(SCFRecipe r, double kT, bool molecule)
{
    if (r!=SCFRecipe::MagneticInsulator && r!=SCFRecipe::Metal) return;
    const std::string name=NameOf(r);
    if (molecule)
        policy::Refuse("SCFRecipe", "scf.recipe "+name+" on a molecule",
                       "a molecule takes Auto or Insulator only: the "+name+" recipe is Kerker-preconditioned Pulay, which mixes in G space "
                       "(doc/SCFRobustnessPlan.md R4.1, ruling (c))");
    if (kT==0.0)
        policy::Refuse("SCFRecipe", "scf.recipe "+name+" at smearingkT 0",
                       "no integer aufbau at a degenerate d frontier or a metal's Fermi level (doc/SCFRobustnessPlan.md §3 row 2': J5, F4); "
                       "state smearingkT >= the frontier splitting (doc/Pins.md pin 15a), typically 1e-3...5e-3 Ha -- smearing is a model "
                       "choice, so the recipe never sets it");
}

//! The accelerator a stage runs, and whether a Ladder must START on its GDM rung (otherwise the caller's accelerator options stand).
struct RecipeAccelerator
{
    SCFAccelerators::Type type;
    bool forceGDMRung = false;   //!< true: a restored kT = 0 Ladder enters at GDM (NaF 9 -> 2 iterations, CLAUDE.md)
};
//! \brief THE ACCELERATOR (§3's accelerator column).  A stated \a asked (not Auto) is obeyed.  Otherwise:
//!   a COLD or CARRIED stage -> Null, so the mixer runs its density history alone (R1.1: no Fock history beside it -- they fight,
//!                              measured 2026-10-22).  A carried stage may have a NEW Hamiltonian or kT: GDM there false-converges on
//!                              its own small steps (F3);
//!   a RESTORED stage at kT = 0 with every virtual built -> the Ladder on its GDM rung (same Hamiltonian, the density IS the answer);
//!   a RESTORED smeared stage   -> DIIS (GDM would minimise an integer D' and never apply the smearing; `SolidState.SmearedExactResume_
//!                                 StartsAtTheSavedState`); a RESTORED partial eigensolve (\a numBands > 0) -> DIIS (GDM needs every virtual).
//! The regime does not enter in B5: every regime's columns name the same accelerators (rows 2 and 3 are always smeared).
inline RecipeAccelerator ChooseAccelerator(SCFRecipe, SCFStageKind k, SCFAccelerators::Type asked, double kT, int numBands)
{
    using T=SCFAccelerators::Type;
    if (asked!=T::Auto)               return {asked};
    if (k!=SCFStageKind::Restored)    return {T::Null};
    if (kT==0.0 && numBands==0)       return {T::Ladder, true};
    return {T::DIIS};
}

//! A stage's SCF parameters with the recipe's mixer knobs filled, and the stated knobs that differ from the recipe (the DEVIATIONS).
struct RecipeKnobs
{
    SCFParams                params;
    std::vector<std::string> deviations;
};
//! \brief THE MIXER KNOBS (§3's mixer column) a NAMED regime fills on a cold or carried stage.  Insulator (and so Auto) fills nothing:
//! its Pulay-D is \c theAutoMixing, which \c ResolveMixing applies to an Auto mixer with no Fock history.  MagneticInsulator and Metal
//! fill \c theKerkerPulayRecipe into every knob left at its \c SCFParams default; a knob the caller stated WINS and is a deviation.  The
//! \c Mixer stays as stated (Auto keeps resolving through \c ResolveMixing, which keeps a hand-set history), so a deck that moves its
//! K+P knobs into `scf.recipe` runs bit-identically.  A RESTORED stage fills nothing (row 4 is the same for every regime).
inline RecipeKnobs FillRecipeKnobs(SCFRecipe r, SCFStageKind k, const SCFParams& p)
{
    RecipeKnobs out{p, {}};
    const SCFRecipe R=Resolved(r);
    if (k==SCFStageKind::Restored || R==SCFRecipe::Insulator) return out;
    const SCFParams d{};
    auto num=[](double x){std::ostringstream s; s<<x; return s.str();};
    if (p.Mixer!=SCFParams::Mixing::Auto)
    {
        out.deviations.push_back("mixer stated (the recipe's knobs are not filled)");
        return out;
    }
    const MixingRecipe& K=theKerkerPulayRecipe;
    SCFParams& q=out.params;
    if (q.KerkerG0==d.KerkerG0)               q.KerkerG0=K.G0;           else if (q.KerkerG0!=K.G0)          out.deviations.push_back("kerkerG0 "+num(q.KerkerG0));
    if (q.PulayDepth==d.PulayDepth)           q.PulayDepth=K.depth;      else if (q.PulayDepth!=K.depth)     out.deviations.push_back("pulayDepth "+num(q.PulayDepth));
    if (q.PulayStart==d.PulayStart)           q.PulayStart=K.start;      else if (q.PulayStart!=K.start)     out.deviations.push_back("pulayStart "+num(q.PulayStart));
    if (q.StartingRelaxRo==d.StartingRelaxRo) q.StartingRelaxRo=K.relax; else if (q.StartingRelaxRo!=K.relax) out.deviations.push_back("startingRelaxRo "+num(q.StartingRelaxRo));
    if (q.UseMOM)                                                        out.deviations.push_back("useMOM true (the recipe runs no MOM: it fights K+P, A7)");
    return out;
}

//! \brief ROW 4's STEP: an Auto accelerator on a RESTORED density pairs with a LINEAR mixer (beside its Fock history), whose START is
//! capped at \c theRestoredLinearStep.  \a resolvesLinear = the stage's mixer resolves to linear (the iterator's \c ResolveMixing).
//! A stated accelerator keeps the caller's step.
inline SCFParams CapRestoredStep(SCFStageKind k, SCFAccelerators::Type asked, const SCFParams& p, bool resolvesLinear)
{
    SCFParams q=p;
    if (k==SCFStageKind::Restored && asked==SCFAccelerators::Type::Auto && resolvesLinear)
        q.StartingRelaxRo=std::min(q.StartingRelaxRo, theRestoredLinearStep);
    return q;
}

//! The §3 row a (regime, stage) cell is.
inline const char* RowOf(SCFRecipe r, SCFStageKind k)
{
    if (k==SCFStageKind::Restored) return "row 4";
    if (k==SCFStageKind::Carried)  return "row 5";
    switch (Resolved(r))
    {
    case SCFRecipe::MagneticInsulator: return "row 2";
    case SCFRecipe::Metal:             return "row 3";
    default:                           return "row 1";
    }
}

//! The mixer a RESOLVED \a p builds, in words.
inline std::string DescribeMixer(const SCFParams& p)
{
    using M=SCFParams::Mixing;
    std::ostringstream s;
    switch (p.Mixer)
    {
    case M::Linear: s<<"linear D-mixing, step "<<p.StartingRelaxRo<<" (ceiling "<<p.MaxRelaxRo<<")"; break;
    case M::Kerker: s<<"Kerker(G0 "<<p.KerkerG0<<"), step "<<p.StartingRelaxRo; break;
    case M::Pulay:
        s<<(p.KerkerG0>0.0 ? "Kerker(G0 " : "Pulay on D (");
        if (p.KerkerG0>0.0) s<<p.KerkerG0<<")-Pulay (";
        s<<"alpha "<<p.StartingRelaxRo<<", depth "<<p.PulayDepth<<", start "<<p.PulayStart<<")"; break;
    case M::Auto:   s<<"Auto (unresolved)"; break;
    }
    return s.str();
}

//! \brief THE POLICY'S REPORT, one per Converge (CLAUDE.md "every policy reports itself"): one console line `[<label> policy] SCFRecipe:
//! ...` and one `policy` row at the caller's report cursor.  \a built = the accelerator \c ChooseAccelerator gave; \a resolved = the
//! stage's parameters as the iterator resolved them; \a statedBy = who stated a non-Auto accelerator ("" = the caller, a deviation).
inline void ReportRecipe(const std::string& label, SCFRecipe r, SCFStageKind k, SCFAccelerators::Type asked, SCFAccelerators::Type built,
                         const SCFParams& resolved, const std::vector<std::string>& deviations, const std::string& statedBy={},
                         const std::string& rowOverride={})
{
    const SCFRecipe R=Resolved(r);
    const std::string acc=NameOf(built);
    std::ostringstream d;
    d<<NameOf(r);
    if (r==SCFRecipe::Auto) d<<" -> "<<NameOf(R);
    d<<", "<<NameOf(k)<<" stage";
    if (k==SCFStageKind::Restored) d<<", kT "<<resolved.SmearingkT;
    const std::string row = rowOverride.empty() ? std::string("SCFRobustnessPlan §3 ")+RowOf(r,k) : rowOverride;
    d<<" ("<<row<<"): accelerator "<<acc;
    if (asked!=SCFAccelerators::Type::Auto) d<<(statedBy.empty() ? " (stated)" : " (stated by "+statedBy+")");
    d<<" + "<<DescribeMixer(resolved);
    if (!deviations.empty())
    {
        d<<"; stated beside it:";
        for (size_t i=0; i<deviations.size(); ++i) d<<(i ? ", " : " ")<<deviations[i];
        d<<" (deviation"<<(deviations.size()>1 ? "s" : "")<<")";
    }
    const std::string why = r==SCFRecipe::Auto
        ? "Auto resolves every stage to Insulator until it can detect a metal or a magnetic frontier (R4.1 as ruled for B5)"
        : "the deck's scf.recipe";
    std::cout<<"["<<label<<" policy] SCFRecipe: "<<d.str()<<"; why: "<<why<<std::endl;
    report::json detail{{"recipe",NameOf(r)}, {"regime",NameOf(R)}, {"stage",NameOf(k)}, {"row",row},
                        {"accelerator",acc}, {"mixer",DescribeMixer(resolved)}, {"deviations",deviations}};
    policy::Report("SCFRecipe", d.str(), policy::Verdict::Choose, why, std::nullopt, detail);
}

} // namespace qchem
