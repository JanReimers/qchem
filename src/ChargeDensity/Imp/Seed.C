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
#include <string>
#include <functional>   // the per-Z valence-count source of the formal-charge rule
#include <optional>     // the run's stated moment (the polarized molecular seed)

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

//! THE FORMAL CHARGE ORDER of an ionic seed, the ONE rule every structure shares (user, 2026-11-12: "we should strive to make seeding
//! ... structure neutral"): \a stated when given, else \c IonicFormalCharges over each site's neutral valence count \a nvalOf(Z).
//! The structure's TYPE never enters; only which library holds its densities does (the periodic seed reads the valence library, the
//! molecular one the all-electron library), and each passes its own \a nvalOf.  \a nval receives the per-site valence counts.
static std::vector<int> SiteFormalCharges(const Structure* st, const std::function<int(int)>& nvalOf, const std::vector<int>& stated,
                                          std::vector<int>& nval)
{
    assert(st);
    std::vector<std::pair<int,int>> atoms;   // {Z, N_val} per atom
    std::map<int,int> nvalByZ;
    st->ForEachSite([&](int Z, const rvec3_t&, bool)
    {
        if (!nvalByZ.count(Z)) nvalByZ[Z]=nvalOf(Z);
        atoms.emplace_back(Z, nvalByZ[Z]);
    });
    nval.clear(); for (const auto& a : atoms) nval.push_back(a.second);
    if (!stated.empty() && stated.size()!=atoms.size())
        throw std::invalid_argument("SiteFormalCharges: "+std::to_string(stated.size())+" stated charges for "+std::to_string(atoms.size())+" sites");
    return stated.empty() ? IonicFormalCharges(atoms) : stated;
}

//! A neutral atom's density AMPLITUDE-SCALED to a charge state the library does not hold (the molecular IonicSAD fallback, the
//! all-electron twin of \c SeedCD's scaled-neutral class).
class ScaledAtomicDensity : public virtual ScalarFunction<double>
{
public:
    ScaledAtomicDensity(std::shared_ptr<const ScalarFunction<double>> f, double scale) : itsF(std::move(f)), itsScale(scale) {}
    virtual double  operator()(const rvec3_t& r) const override {return itsScale*(*itsF)(r);}
    virtual rvec3_t Gradient  (const rvec3_t& r) const override {return itsScale*itsF->Gradient(r);}
private:
    std::shared_ptr<const ScalarFunction<double>> itsF;
    double itsScale;
};

//! The neutral valence electron count of \a Z from the periodic table's valence configuration (the all-electron library's twin of
//! the valence library's neutral charge).
static int PeriodicValence(int Z)
{
    const size_t* vc=thePeriodicTable().GetValanceConfiguration(size_t(Z));
    return int(vc[0]+vc[1]+vc[2]+vc[3]);
}

//! The plane-wave seed's own density-fit basis, from the WHOLE-SET fit factory.  It is not the Hamiltonian's: its grid reports one
//! level into a context naming the seed.
static std::shared_ptr<const BasisSet::cFIT_CD_ABS> SeedFit(const BasisSet::tBasisSet<dcmplx>* bs, const Structure* st)
{
    qchem::report::Scope here("seedDensity");
    return std::shared_ptr<const BasisSet::cFIT_CD_ABS>(bs->CreateCDFitBasisSet(st));
}

//! rho(r) = N/V on the first block: the seed both lineages share.
template <class T> static tChargeDensity<T>* UniformSeed(const BasisSet::tBasisSet<T>* bs, const ElectronConfiguration* ec, const DensityMaker& maker)
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
            comp->Insert(maker.Make<double>(D0, rb, irr), irr);
            return comp;
        }
    const tobs_t<T>* block = (*bs)[0];
    const size_t     n     = block->GetNumFunctions();
    const int        N     = ec->GetN(irr);
    hmat_t<T> D0=blazem::zeroH<T>(n);
    for (size_t i=0;i<n;i++) D0(i,i)=T(double(N)/double(n));
    return maker.Make<T>(D0, block, irr).release();
}

//! One site of the molecular SAD superposition: its centre, its total radial density at its electron count with the amplitude scale
//! (1, or the scaled-neutral fallback's target/base), and its Hund pair at that count when the library holds one.
struct SADSite
{
    rvec3_t R;
    std::shared_ptr<const RadialDensity> rho;            //!< the site's total
    double scale=1.0;                                    //!< amplitude scale on rho (and on the pair)
    std::shared_ptr<const RadialDensity> maj, min;       //!< the Hund pair (majority, minority), null when the library holds none
    int spin=0;                                          //!< +1 / -1 (a flipped site) / 0 (no planted moment)
};

//! The site's total as the unpolarized seed has always inserted it: the recentred table, amplitude-scaled only when the scale is not 1.
static std::shared_ptr<const ScalarFunction<double>> SiteTotal(const SADSite& site)
{
    auto rc=std::make_shared<const RecentredAtomicDensity>(site.rho, site.R);
    if (site.scale==1.0) return rc;
    return std::make_shared<const ScaledAtomicDensity>(rc, site.scale);
}

//! The spin-agnostic seed: the sum of every site's total (the molecular SAD seed as it was before the polarized one).
static NumericCD* TotalSeed(const std::vector<SADSite>& sites, double N, std::shared_ptr<const qcMesh::Mesh> mesh)
{
    auto* cd=new NumericCD(N, std::move(mesh));
    for (const auto& site : sites) cd->Insert(SiteTotal(site));
    return cd;
}

//! The run's net moment nUp - nDown when its configuration FIXES each channel's count; nullopt when the channels share one reservoir
//! (the moment is an output) or the basis has no spin channels.
static std::optional<double> StatedMoment(const BasisSet::tBasisSet<double>* bs, const ElectronConfiguration* ec)
{
    const ReservoirPartition part=ec->GetPartition();
    if (part.spansSpin) return std::nullopt;
    const auto count=[&](Spin s)
    {
        const auto irs=bs->GetIrreps(s);
        if (irs.empty()) return -1;
        if (part.spansSpatial) return ec->GetN(irs[0]);
        int n=0; for (const auto& irr : irs) n+=ec->GetN(irr);
        return n;
    };
    const int up=count(Spin::Up), dn=count(Spin::Down);
    if (up<0 || dn<0) return std::nullopt;
    return double(up-dn);
}

//! THE POLARIZED molecular seed (see \c MakeSeedDensity): each channel the sum of the sites' channel densities.  A site's magnetization
//! m = maj - min (sign = its spin) is scaled by ONE factor f so the planted net moment is the run's \a moment (when stated and the plant
//! is not zero); its channels are rho/2 +- f m/2 = a maj + b min with a = (1+f)/2, b = (1-f)/2, swapped on a flipped site.  A pairless
//! site splits rho/2.  Returns the spin-agnostic seed when nothing is planted: no pair, or a ferromagnetic plant scaled to a stated 0.
static tChargeDensity<double>* PolarizedSeed(const std::vector<SADSite>& sites, double N, std::optional<double> moment,
                                             std::shared_ptr<const qcMesh::Mesh> mesh)
{
    double planted=0; bool anyPair=false;
    for (const auto& site : sites)
        if (site.maj && site.spin!=0) { anyPair=true; planted+=site.spin*site.scale*(site.maj->Charge()-site.min->Charge()); }
    // A planted 0 (an AFM pattern of flipped sites, a broken-symmetry singlet) keeps its pattern: no factor maps 0 to the stated moment.
    const double f=(moment && std::abs(planted)>1e-8) ? *moment/planted : 1.0;
    if (!anyPair || f==0.0) return TotalSeed(sites, N, std::move(mesh));   // nothing to plant: the spin-agnostic seed, as before
    const double a=0.5*(1+f), b=0.5*(1-f);
    // The channels' counts (N +- M)/2 take the STATED moment itself, not f*planted (= it to an ulp): a one-electron doublet's minority
    // channel is then exactly 0, never -1e-17 (D-SEED-EMPTY-CHANNEL).
    const double M=(moment && std::abs(planted)>1e-8) ? *moment : planted;   // f*planted, exactly
    auto up=std::make_unique<NumericCD>(0.5*(N+M), mesh);
    auto dn=std::make_unique<NumericCD>(0.5*(N-M), mesh);
    const auto scaled=[](std::shared_ptr<const RadialDensity> r, const rvec3_t& R, double c) -> std::shared_ptr<const ScalarFunction<double>>
    {
        return std::make_shared<const ScaledAtomicDensity>(std::make_shared<const RecentredAtomicDensity>(std::move(r), R), c);
    };
    for (const auto& site : sites)
    {
        if (!site.maj || site.spin==0)
        {
            up->Insert(scaled(site.rho, site.R, 0.5*site.scale));
            dn->Insert(scaled(site.rho, site.R, 0.5*site.scale));
            continue;
        }
        const auto& hi=site.spin>0 ? site.maj : site.min;   // the channel this site's majority lands in is Up unless it is flipped
        const auto& lo=site.spin>0 ? site.min : site.maj;
        up->Insert(scaled(hi, site.R, a*site.scale)); up->Insert(scaled(lo, site.R, b*site.scale));
        dn->Insert(scaled(lo, site.R, a*site.scale)); dn->Insert(scaled(hi, site.R, b*site.scale));
    }
    return new PolarizedNumericCD(std::move(up), std::move(dn));
}

tChargeDensity<double>* MakeSeedDensity(SeedStrategy s, const BasisSet::tBasisSet<double>* bs,
                                        const Structure* st, const ElectronConfiguration* ec, const qcMesh::MeshParams& quadrature,
                                        const qcMesh::MeshMaker& meshMaker, const SeedKey& key, const DensityMaker& maker, bool polarized)
{
    assert(bs);
    assert(ec);
    if (s==SeedStrategy::Default || s==SeedStrategy::CoreGuess) return nullptr;   // density-independent core-Hamiltonian guess (cd=0)
    if (s==SeedStrategy::Uniform) return UniformSeed<double>(bs, ec, maker);
    if (s!=SeedStrategy::SAD && s!=SeedStrategy::IonicSAD)
        throw std::runtime_error("MakeSeedDensity: unknown SeedStrategy");   // see the IonicSAD arm (V1.19)
    if (!st) throw std::logic_error("MakeSeedDensity: the SAD/IonicSAD seed needs a Structure (atom Z + positions)");

    // SAD: the superposition of neutral atomic densities; each atom is the RUN's atom (SeedLibrary): the pseudo-atom at its q under a
    // pseudopotential, the all-electron atom otherwise.  IonicSAD: the same at each site's FORMAL CHARGE (Li+ F-), from the same charge
    // rule as the periodic seed (SiteFormalCharges): a site's electron count is its pseudo-atom's (the valence density's charge at the
    // key's q) under a pseudopotential, Z all-electron; a charge state the library lacks is the neutral density amplitude-scaled to it
    // (SeedCD's fallback); a stripped site (H+) carries nothing.  With every charge zero IonicSAD IS the SAD seed.  The molecule's own net
    // charge is not distributed: the fit's charge constraint carries it.
    const auto neutral=[&key](int Z){ return GetAtomicDensity(Z, key.functional, SeedLibrary(key,Z), -1, key.Q(Z)); };
    std::vector<int> q;   // per-site formal charges (IonicSAD only)
    if (s==SeedStrategy::IonicSAD)
    {
        // The electronegativity rule moves VALENCE electrons, whichever library: a pseudo-atom's whole count, an all-electron atom's
        // periodic-table valence.
        const auto nvalOf=[&](int Z){ return key.Q(Z)>=0 ? int(std::lround(neutral(Z).Charge())) : PeriodicValence(Z); };
        std::vector<int> nval;
        q=SiteFormalCharges(st, nvalOf, {}, nval);
    }
    std::vector<SADSite> sites;
    size_t i=0;
    st->ForEachSite([&](int Z, const rvec3_t& R, bool flip)
    {
        const size_t site=i++;
        const std::string lib=SeedLibrary(key,Z);
        SADSite one{.R=R};
        int target=-1;   // the site's electron count; -1 = the neutral (the library's first match)
        if (s==SeedStrategy::IonicSAD)
        {
            // The site's own count (the base its charge is taken from) is the pseudo-atom's, or Z.
            const int base=key.Q(Z)>=0 ? int(std::lround(neutral(Z).Charge())) : Z;
            target=base-q[site];
            if (target<=0) return;
            if (!HasAtomicDensity(Z, key.functional, lib, target, key.Q(Z)))
            {
                one.rho=std::make_shared<const RadialDensity>(neutral(Z));
                one.scale=double(target)/double(base);
                target=-1;   // the scaled NEUTRAL: its pair (if any) is the neutral's, scaled alike
            }
        }
        if (!one.rho) one.rho=std::make_shared<const RadialDensity>(GetAtomicDensity(Z, key.functional, lib, target, key.Q(Z)));
        if (polarized && HasAtomicSpinPair(Z, key.functional, lib, target, key.Q(Z)))
        {
            auto [maj, min]=GetAtomicSpinPair(Z, key.functional, lib, target, key.Q(Z));
            one.maj=std::make_shared<const RadialDensity>(std::move(maj));
            one.min=std::make_shared<const RadialDensity>(std::move(min));
            one.spin=key.siteSpins.empty() ? (flip ? -1 : +1) : key.siteSpins.at(site);   // a STATED ordering wins over the flips
        }
        sites.push_back(std::move(one));
    });
    // The seed is only a function on space, so it OWNS the points it integrates itself on (a Coulomb fit has none, doc/Pins.md):
    // the run's XC mesh, built by the run point's maker and reported one level into the seed's own context.
    auto mesh=[&]{ qchem::report::Scope here("seedDensity"); return std::make_shared<const qcMesh::Mesh>(meshMaker.Make(quadrature)); }();
    const double N=st->GetNumElectrons();
    if (!polarized) return TotalSeed(sites, N, std::move(mesh));
    return PolarizedSeed(sites, N, StatedMoment(bs, ec), std::move(mesh));
}

tChargeDensity<dcmplx>* MakeSeedDensity(SeedStrategy s, const BasisSet::tBasisSet<dcmplx>* bs,
                                        const UnitCell& cell, const ElectronConfiguration* ec, const DensityMaker& maker,
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
        return UniformSeed<dcmplx>(bs, ec, maker);
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
    // A STATED charge order (the run's siteCharges, checked at resolve) wins; else the heuristic: Na+1, F-1, charge conserved.
    std::vector<int> nval;
    const std::vector<int> q=SiteFormalCharges(st,
        [&key](int Z){ return (int)std::lround(GetAtomicDensity(Z,key.functional,"atomic_valence_densities.json",-1,key.Q(Z)).Charge()); },
        key.siteCharges, nval);
    std::vector<int> targets(nval.size());                    // site i -> TARGET valence count N_val-q_i
    for (size_t i=0;i<nval.size();i++)
        targets[i] = nval[i] - q[i];
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
