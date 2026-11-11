// File: ChargeDensity/Imp/Seed.C  SCF seed-density factory implementation (Phase 0: CoreGuess / Uniform).
module;
#include <cassert>
#include <type_traits>
#include <cstddef>
#include <memory>
#include <stdexcept>   // the unsupported (strategy x T) arms THROW (V1.19: an assert-only arm is a silent core guess in Release)
#include <vector>
#include <utility>
#include <algorithm>
#include <map>
#include <cmath>

module qchem.ChargeDensity.Seed;
import qchem.PeriodicTable;                   // thePeriodicTable -> GetElectronegativity (the ionic heuristic)
import qchem.ChargeDensity.Factory;          // IrrepCD_Factory<T>
import qchem.ChargeDensity.Types;            // tobs_t<T>
import qchem.CompositeCD;                    // tComposite_CD<dcmplx> (the mixed Uniform seed's run-shaped wrap)
import qchem.ChargeDensity.NumericCD;// NumericCD (the molecular SAD seed, double only)
import qchem.ChargeDensity.SeedCD;    // SeedCD (the plane-wave SAD seed, dcmplx only)
import qchem.ChargeDensity.AtomicDensity;    // GetAtomicDensity, RadialDensity, RecentredAtomicDensity
import qchem.BasisSet.Orbital_DFT_IBS;           // Orbital_DFT_IBS<dcmplx> (the PW block: CreateCDFitBasisSet for the seed)
import qchem.Structure;                       // Structure, Atom (atom Z + positions)
import qchem.Blaze;                           // blazem::zeroH, hmat_t
import qchem.Types;                           // dcmplx
import qchem.Reporting;                       // report::Scope (SeedFit)

namespace qchem::ChargeDensity
{

std::vector<int> IonicFormalCharges(const std::vector<std::pair<int,int>>& atoms)
{
    const size_t n=atoms.size();
    std::vector<int>    q(n, 0);
    std::vector<double> en(n);
    std::vector<int>    giveCap(n), takeCap(n);   // electrons this atom can still donate / accept
    const int closed[]={2,8,18,32};               // closed-shell valence-electron counts an acceptor fills to
    for (size_t i=0;i<n;i++)
    {
        int Z=atoms[i].first, nval=atoms[i].second;
        en[i]      = thePeriodicTable().GetElectronegativity(Z);
        giveCap[i] = nval;                          // a donor can shed its whole valence (down to the noble core)
        int ceil=nval; for (int c:closed) if (c>=nval) { ceil=c; break; }
        takeCap[i] = ceil-nval;                     // an acceptor fills up to the next closed shell
    }
    // Greedy charge-conserving transfer: move one electron at a time from the lowest-EN atom that can still
    // donate to the highest-EN atom that can still accept, while a real EN gap remains.
    // TIES (equal EN = the same element) go to the atom with the most capacity LEFT, so equivalent atoms share
    // the transfer evenly.  A strict EN comparison alone drained the FIRST donor's whole valence before touching
    // the next: beta-MnO2 came out Mn+7 / Mn+1 instead of Mn4+ / Mn4+ (2026-10-22), seeding one Mn as nothing and
    // the other non-magnetic, so the A1-AFM order never formed.
    const double kGap=0.5;                          // min Pauling gap to call a pair ionic (NaF 3.0, CsI 1.9)
    const double kTie=1e-9;
    while (true)
    {
        int d=-1, a=-1;
        for (size_t i=0;i<n;i++)
        {
            if (giveCap[i]>0 && (d<0 || en[i]<en[d]-kTie || (std::abs(en[i]-en[d])<=kTie && giveCap[i]>giveCap[d]))) d=(int)i;
            if (takeCap[i]>0 && (a<0 || en[i]>en[a]+kTie || (std::abs(en[i]-en[a])<=kTie && takeCap[i]>takeCap[a]))) a=(int)i;
        }
        if (d<0 || a<0 || d==a || en[a]-en[d]<kGap) break;
        q[d]+=1; q[a]-=1; giveCap[d]-=1; takeCap[a]-=1;   // donor +charge, acceptor -charge; charge conserved
    }
    return q;
}

//! The plane-wave seed's own density-fit basis, from the WHOLE-SET fit factory.  It is not the Hamiltonian's: its grid reports one
//! level into a context naming the seed.
static std::shared_ptr<const BasisSet::cFIT_CD_ABS> SeedFit(const BasisSet::tBasisSet<dcmplx>* bs, const Structure* st)
{
    qchem::report::Scope here("seedDensity");
    return std::shared_ptr<const BasisSet::cFIT_CD_ABS>(bs->CreateCDFitBasisSet(st));
}

//! rho(r) = N/V on the first block: the seed both lineages share.
template <class T> static tChargeDensity<T>* UniformSeed(const BasisSet::tBasisSet<T>* bs, const ElectronConfiguration* ec)
{
    // rho(r) = N/V  =>  D = (N/n) I on the first block.  A single block suffices since rho is
    // constant: every block's first Hartree/XC build sees the same total density.  This is the
    // centralized version of the old per-test plane-wave seed boilerplate.
    // MIXED-AWARE (doc/RealComplexPlan.md 3c-3): after the factory flip a REAL TRIM block may sit
    // at index 0 (Γ-first), where the same-scalar [0] view throws -- its seed is the <double> leaf
    // inside the run-shaped composite (the same uniform rho either way).
    const Irrep irr = bs->GetIrreps(Spin::None)[0];
    if constexpr (std::is_same_v<T,dcmplx>)
        if (const auto* rb = bs->GetRealIBS(0))
        {
            const size_t n = rb->GetNumFunctions();
            const int    N = ec->GetN(irr);
            hmat_t<double> D0=blazem::zeroH<double>(n);
            for (size_t i=0;i<n;i++) D0(i,i)=double(N)/double(n);
            // The run's density group, like every composite of this run (SolidState): the seed's ρ̃ AND its raster are star-averaged
            // under it.  Before 2026-11-09 this composite held no ops, yet its raster was averaged anyway (the fit basis carried
            // them) while its ρ̃ was not -- one density, two symmetries (V-SYM-SEQUESTER-OPS).
            const auto* imp=bs->GetImposition();
            auto* comp=new tComposite_CD<dcmplx>(imp ? imp->DensityGroup() : nullptr);
            comp->Insert(std::unique_ptr<tDM_CD<double>>(IrrepCD_Factory<double>(D0, rb, irr)), irr);
            return comp;
        }
    const tobs_t<T>* block = (*bs)[0];
    const size_t     n     = block->GetNumFunctions();
    const int        N     = ec->GetN(irr);
    hmat_t<T> D0=blazem::zeroH<T>(n);
    for (size_t i=0;i<n;i++) D0(i,i)=T(double(N)/double(n));
    return IrrepCD_Factory<T>(D0, block, irr);
}

tChargeDensity<double>* MakeSeedDensity(SeedStrategy s, const BasisSet::tBasisSet<double>* bs,
                                        const Structure* st, const ElectronConfiguration* ec, const qcMesh::MeshParams& quadrature)
{
    assert(bs);
    assert(ec);
    switch (s)
    {
    case SeedStrategy::Default:   // the molecular present-day behaviour
    case SeedStrategy::CoreGuess:
        return nullptr;   // density-independent core-Hamiltonian guess (today's cd=0)
    case SeedStrategy::Uniform:
        return UniformSeed<double>(bs, ec);
    case SeedStrategy::SAD:
    {
        // Superposition of neutral atomic densities.  Drop each atom's recentred radial density (LDA DB) into a NumericCD;
        // FittedVee/FittedVxc consume it (Coulomb projection + rho(r)).
        if (!st) throw std::logic_error("MakeSeedDensity: the SAD seed needs a Structure (atom Z + positions)");
        // The seed is only a function on space, so it OWNS the points it integrates itself on (a Coulomb fit has none, doc/Pins.md):
        // the run's XC mesh, reported one level into the seed's own context.
        auto mesh=[&]{ qchem::report::Scope here("seedDensity"); return std::make_shared<const qcMesh::Mesh>(st->CreateIntegrationMesh(quadrature)); }();
        auto* cd = new NumericCD(st->GetNumElectrons(), std::move(mesh));
        st->ForEachSite([&cd](int Z, const rvec3_t& R, bool)
        {
            auto rad = std::make_shared<const RadialDensity>(GetAtomicDensity(Z));
            cd->Insert(std::make_shared<const RecentredAtomicDensity>(rad, R));
        });
        return cd;
    }
    case SeedStrategy::IonicSAD:
        // THROW, not assert: an assert-only arm returns null under -DNDEBUG and SCFIterator silently
        // runs a core guess instead of the seed the caller asked for (V1.19).
        throw std::runtime_error("MakeSeedDensity: molecular (AO) IonicSAD is a later phase; "
                                 "use SAD or the plane-wave path");
    default:
        throw std::runtime_error("MakeSeedDensity: unknown SeedStrategy");   // see the IonicSAD arm (V1.19)
    }
}

tChargeDensity<dcmplx>* MakeSeedDensity(SeedStrategy s, const BasisSet::tBasisSet<dcmplx>* bs,
                                        const UnitCell& cell, const ElectronConfiguration* ec,
                                        bool polarized, const SeedKey& key)
{
    assert(bs);
    assert(ec);
    switch (s)
    {
    case SeedStrategy::CoreGuess:
        return nullptr;   // density-independent core-Hamiltonian guess (today's cd=0)
    case SeedStrategy::Default:   // the plane-wave present-day behaviour
    case SeedStrategy::Uniform:
        return UniformSeed<dcmplx>(bs, ec);
    case SeedStrategy::SAD:
    {
        // Plane-wave (FT) SAD: rho-tilde(G) = Sum_atoms rho_atom(|G|) e^{-iG.R}, assembled by the basis.
        // A polarized run gets the two-channel seed (§10): Hund pairs + per-atom flips choose the basin.
        // The WHOLE-SET fit factory (mixed-aware since 3c-3) -- serves from the first block of
        // either scalar, so a Γ-first real TRIM block no longer needs a block-0 cast here.
        auto fb=SeedFit(bs, &cell);
        if (polarized) return new PolarizedSeedCD(fb, cell, key);
        return new SeedCD(fb, cell, key);
    }
    case SeedStrategy::IonicSAD:
    {
        // Ionic superposition: like SAD but each species' valence density is scaled to its formal-charge
        // electron count (Na+ -> 0, F- -> 8/7), pre-baking the charge transfer the neutral SAD leaves for
        // the SCF.  The formal charges come from the structure-wide electronegativity heuristic
        // (IonicFormalCharges); per-species N_val is the valence density's charge.
        // Whole-set fit factory, as in the SAD arm above (mixed-aware since 3c-3).
        auto fb=SeedFit(bs, &cell);

        const std::vector<int> targets = IonicSADTargets(&cell, key);   // the ONE resolution (S3 shares it), per site
        if (polarized) return new PolarizedSeedCD(fb, cell, key, targets);   // §10: the two-channel ionic seed
        return new SeedCD(fb, cell, key, targets);               // SeedCD prefers the DIFFUSE charge-state density
    }
    default:
        throw std::runtime_error("MakeSeedDensity: unknown SeedStrategy");   // see the IonicSAD arm (V1.19)
    }
}

// The IonicSAD per-SITE TARGET valence counts -- the ONE resolution the seed AND the S3 magnetic
// decoration share (formal charges from the electronegativity heuristic; N_val from the neutral
// valence density's charge): F -> 7-(-1)=8 (F-), Na -> 1-1=0 (Na+).  One entry per atom: a per-species map
// kept only the LAST atom's target, collapsing a mixed-valence cell onto one charge state.
std::vector<int> IonicSADTargets(const Structure* st, const SeedKey& key)
{
    assert(st);
    std::vector<std::pair<int,int>> atoms;                      // {Z, N_val} per atom
    std::map<size_t,int> nvalByZ;
    st->ForEachSite([&](int Z, const rvec3_t&, bool)
    {
        if (!nvalByZ.count(Z))
            nvalByZ[Z] = (int)std::lround(GetAtomicDensity(Z,key.functional,"atomic_valence_densities.json",-1,key.Q(Z)).Charge());
        atoms.emplace_back(Z, nvalByZ[Z]);
    });
    // A STATED charge order (the run's siteCharges, checked at resolve) wins; else the heuristic: Na+1, F-1, charge conserved.
    assert(key.siteCharges.empty() || key.siteCharges.size()==atoms.size());
    const std::vector<int> q = key.siteCharges.empty() ? IonicFormalCharges(atoms) : key.siteCharges;
    std::vector<int> targets(atoms.size());                    // site i -> TARGET valence count N_val-q_i
    for (size_t i=0;i<atoms.size();i++)
        targets[i] = atoms[i].second - q[i];
    return targets;
}

// Shubnikov S3: the per-atom collinear labels, resolved by the SEED's own rule (SeedCD's ctor), SITE by site:
// target = the site's ionic entry, else the neutral valence count; magnetic iff a spin pair exists AT that
// target (a stripped cation, target<=0, is non-magnetic).  Label = the atom's flip bit.
std::vector<int> MagneticDecoration(const Structure* st, const SeedKey& seedKey,
                                    const std::vector<int>& ionicNvalBySite)
{
    assert(st);
    assert(ionicNvalBySite.empty() || ionicNvalBySite.size()==st->GetNumAtoms());
    const std::string db = "atomic_valence_densities.json";
    std::map<std::pair<int,int>,bool> magnetic;   // (Z, target) -> has a spin pair
    std::vector<int> spins;
    st->ForEachSite([&](int iZ, const rvec3_t&, bool spinFlip)
    {
        const size_t i = spins.size();
        const int target = ionicNvalBySite.empty() ? (int)std::lround(GetAtomicDensity(iZ, seedKey.functional, db, -1, seedKey.Q(iZ)).Charge())
                                                   : ionicNvalBySite[i];
        const auto key = std::make_pair(iZ, target);
        if (!magnetic.count(key)) magnetic[key] = target>0 && HasAtomicSpinPair(iZ, seedKey.functional, db, target, seedKey.Q(iZ));
        spins.push_back(magnetic.at(key) ? (spinFlip ? -1 : +1) : 0);
    });
    return spins;
}

} //namespace
