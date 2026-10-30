// File: Calculation/Imp/MoleculeSetup.C  The molecular run setup: every decision made ONCE, before a geometry point is built.
module;
#include <memory>
#include <stdexcept>
#include <string>
module qchem.MoleculeSetup;

import qchem.PeriodicTable;   // thePeriodicTable (Z -> element symbol, for the PP lookup)

namespace qchem
{

std::shared_ptr<const MoleculeSetup> MoleculeSetup::Resolve(const Structure& st, const CalcOptions& opts, const AcceleratorOptions& acc)
{
    auto s=std::make_shared<MoleculeSetup>();
    s->opts=opts;
    s->acc =acc;
    // An open shell (2S>0) is unrestricted: it needs distinct up/down densities, so promote to polarized.
    if (s->opts.multiplicity > 1) s->opts.spin = SpinGroup::Polarized;
    if (!opts.pp.empty()) s->pp=Pseudopotential::PseudopotentialSet::Get(opts.pp).With(opts.ppOverride);
    s->reference=st.Clone();   // polymorphic deep copy: a UnitCell stays periodic, not sliced to a Molecule
    const std::shared_ptr<const Structure> point=s->PointStructure(*s->reference);
    s->Ne=int(point->GetNumElectrons());

    // The multiplicity 2S+1 -> (nUp, nDown) must be realisable: 0 <= 2S <= Ne and Ne-2S even (fail loud, not a silent miscount).
    if (opts.multiplicity > 0)
    {
        const int twoS = opts.multiplicity - 1;
        if (twoS > s->Ne || ((s->Ne - twoS) % 2) != 0)
            throw std::runtime_error("qchem::Calculation: multiplicity " + std::to_string(opts.multiplicity) +
                " (2S=" + std::to_string(twoS) + ") is incompatible with " + std::to_string(s->Ne) +
                " electrons -- need 0 <= 2S <= Ne and Ne-2S even.");
    }

    // Seed: an explicit opts.seed wins; otherwise auto -- DFT from superposition-of-atomic-densities (SAD), HF/1-e from the
    // core guess (Default).
    using qchem::ChargeDensity::SeedStrategy;
    const bool dft = qchem::Hamiltonian::IsDFT(opts.model);
    s->seed = (opts.seed != SeedStrategy::Default) ? opts.seed : (dft ? SeedStrategy::SAD : SeedStrategy::Default);

    using AType = qchem::SCFAccelerators::Type;
    s->accelerator = acc.type=="DIIS"   ? AType::DIIS
                   : acc.type=="GDM"    ? AType::GDM
                   : acc.type=="Ladder" ? AType::Ladder
                   : throw std::runtime_error("qchem::Calculation: unknown accelerator type \"" + acc.type +
                                              "\" (expected DIIS | GDM | Ladder)");
    // DFT (and the LSDA pseudopotential) need DIIS engaged from iteration 0: the default Z-scaled EMax gate keeps DIIS off and
    // the SCF limit-cycles to a non-converged snapshot (the "M_Sym layout UB" mechanism, see M_DFT).  So default them to a high
    // EMax ("DIIS from start") unless the caller pinned one.
    const bool   dftLike = dft || !opts.pp.empty();
    const int    Z       = point->GetNuclearCharge();
    s->eMax = acc.eMax > 0.0 ? acc.eMax : (dftLike ? 100.0 : Z*Z*0.08/32);
    return s;
}

// Re-express a structure as its VALENCE ions for a pseudopotential run: each atom keeps its true species Z (so the PP lookup and
// the true geometry are intact) but carries net charge Z-Zion(its own element, from the run's chosen set), so GetNumElectrons()
// reports the total Zion valence count the EC + FittedVee charge constraint consume.  Clone-and-mutate (not rebuild-as-Molecule)
// so the concrete geometry is preserved: a PP run on a periodic UnitCell stays a UnitCell (Ewald ion-ion, uniform mesh).
std::shared_ptr<const Structure> MoleculeSetup::PointStructure(const Structure& st) const
{
    auto point = st.Clone();
    if (pp)
        for (size_t a=0; a<point->GetNumAtoms(); a++)
        {
            Atom* atom = (*point)[a];
            atom->itsCharge = double(atom->itsZ - pp->Valence(thePeriodicTable().GetSymbol(atom->itsZ)));   // net charge = Z - Zion
        }
    return point;
}

} // namespace qchem
