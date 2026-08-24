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
import qchem.Mesh;                            // qcMesh::MeshParams (CreateCDFitBasisSet arg)
import qchem.Structure;                       // Structure, Atom (atom Z + positions)
import qchem.Blaze;                           // blazem::zeroH, hmat_t
import qchem.Types;                           // dcmplx

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
    const double kGap=0.5;                          // min Pauling gap to call a pair ionic (NaF 3.0, CsI 1.9)
    while (true)
    {
        int d=-1, a=-1;
        for (size_t i=0;i<n;i++)
        {
            if (giveCap[i]>0 && (d<0 || en[i]<en[d])) d=(int)i;
            if (takeCap[i]>0 && (a<0 || en[i]>en[a])) a=(int)i;
        }
        if (d<0 || a<0 || d==a || en[a]-en[d]<kGap) break;
        q[d]+=1; q[a]-=1; giveCap[d]-=1; takeCap[a]-=1;   // donor +charge, acceptor -charge; charge conserved
    }
    return q;
}

template <class T> tChargeDensity<T>* MakeSeedDensity(SeedStrategy s, const BasisSet::tBasisSet<T>* bs,
                                                      const Structure* st, const ElectronConfiguration* ec,
                                                      bool polarized)
{
    assert(bs);
    assert(ec);
    // Default resolves to each path's present-day behaviour: molecular core guess, plane-wave uniform.
    if (s==SeedStrategy::Default)
        s = std::is_same_v<T,double> ? SeedStrategy::CoreGuess : SeedStrategy::Uniform;

    switch (s)
    {
    case SeedStrategy::CoreGuess:
        return nullptr;   // density-independent core-Hamiltonian guess (today's cd=0)

    case SeedStrategy::Uniform:
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
                auto* comp=new tComposite_CD<dcmplx>();
                comp->Insert(std::unique_ptr<tDM_CD<double>>(IrrepCD_Factory<double>(D0, rb, irr)));
                return comp;
            }
        const tobs_t<T>* block = (*bs)[0];
        const size_t     n     = block->GetNumFunctions();
        const int        N     = ec->GetN(irr);
        hmat_t<T> D0=blazem::zeroH<T>(n);
        for (size_t i=0;i<n;i++) D0(i,i)=T(double(N)/double(n));
        return IrrepCD_Factory<T>(D0, block, irr);
    }

    case SeedStrategy::SAD:
    {
        // Superposition of neutral atomic densities (the molecular SAD seed, double/AO path).  Drop each
        // atom's recentred radial density (LDA DB) into a NumericCD; FittedVee/FittedVxc consume it
        // (Coulomb projection + rho(r)).  The plane-wave (dcmplx, FT) face is Phase 2.
        if constexpr (std::is_same_v<T,double>)
        {
            assert(st && "SAD seed needs a Structure (atom Z + positions)");
            auto* cd = new NumericCD(st->GetNumElectrons());
            st->ForEachSite([&cd](int Z, const rvec3_t& R, bool)
            {
                auto rad = std::make_shared<const RadialDensity>(GetAtomicDensity(Z));
                cd->Insert(std::make_shared<const RecentredAtomicDensity>(rad, R));
            });
            return cd;
        }
        else
        {
            // Plane-wave (FT) SAD: rho-tilde(G) = Sum_atoms rho_atom(|G|) e^{-iG.R}, assembled by the basis.
            // A polarized run gets the two-channel seed (§10): Hund pairs + per-atom flips choose the basin.
            assert(st && "SAD plane-wave seed needs a Structure");
            // The WHOLE-SET fit factory (mixed-aware since 3c-3) -- serves from the first block of
            // either scalar, so a Γ-first real TRIM block no longer needs a block-0 cast here.
            std::shared_ptr<const BasisSet::cFIT_CD_ABS> fb(bs->CreateCDFitBasisSet(st, qcMesh::MeshParams{}));
            if (polarized) return new PolarizedSeedCD(fb, st);
            return new SeedCD(fb, st);
        }
    }
    case SeedStrategy::IonicSAD:
    {
        // Ionic superposition: like SAD but each species' valence density is scaled to its formal-charge
        // electron count (Na+ -> 0, F- -> 8/7), pre-baking the charge transfer the neutral SAD leaves for
        // the SCF.  The formal charges come from the structure-wide electronegativity heuristic
        // (IonicFormalCharges); per-species N_val is the valence density's charge.  Plane-wave path only for
        // now (NaF/CsI are PW); molecular IonicSAD is a later phase.
        if constexpr (std::is_same_v<T,dcmplx>)
        {
            assert(st && "IonicSAD plane-wave seed needs a Structure");
            // Whole-set fit factory, as in the SAD arm above (mixed-aware since 3c-3).
            std::shared_ptr<const BasisSet::cFIT_CD_ABS> fb(bs->CreateCDFitBasisSet(st, qcMesh::MeshParams{}));

            std::map<size_t,int> targetByZ = IonicSADTargets(st, "LDA");   // the ONE resolution (S3 shares it)
            if (polarized) return new PolarizedSeedCD(fb, st, "LDA", targetByZ);   // §10: the two-channel ionic seed
            return new SeedCD(fb, st, "LDA", targetByZ);               // SeedCD prefers the DIFFUSE charge-state density
        }
        else
        {
            // THROW, not assert: an assert-only arm returns null under -DNDEBUG and SCFIterator silently
            // runs a core guess instead of the seed the caller asked for (V1.19).
            throw std::runtime_error("MakeSeedDensity: molecular (AO) IonicSAD is a later phase; "
                                     "use SAD or the plane-wave path");
        }
    }

    default:
        throw std::runtime_error("MakeSeedDensity: unknown SeedStrategy");   // see the IonicSAD arm (V1.19)
    }
}

template tChargeDensity<double>* MakeSeedDensity<double>(SeedStrategy, const BasisSet::tBasisSet<double>*,
                                                        const Structure*, const ElectronConfiguration*, bool);
template tChargeDensity<dcmplx>* MakeSeedDensity<dcmplx>(SeedStrategy, const BasisSet::tBasisSet<dcmplx>*,
                                                        const Structure*, const ElectronConfiguration*, bool);

// The IonicSAD per-species TARGET valence counts -- the ONE resolution the seed AND the S3 magnetic
// decoration share (formal charges from the electronegativity heuristic; N_val from the neutral
// valence density's charge): F -> 7-(-1)=8 (F-), Na -> 1-1=0 (Na+).
std::map<size_t,int> IonicSADTargets(const Structure* st, const std::string& functional)
{
    assert(st);
    std::vector<std::pair<int,int>> atoms;                      // {Z, N_val} per atom
    std::map<size_t,int> nvalByZ;
    st->ForEachSite([&](int Z, const rvec3_t&, bool)
    {
        if (!nvalByZ.count(Z))
            nvalByZ[Z] = (int)std::lround(GetAtomicDensity(Z,functional,"atomic_valence_densities.json").Charge());
        atoms.emplace_back(Z, nvalByZ[Z]);
    });
    std::vector<int> q = IonicFormalCharges(atoms);             // Na+1, F-1; conserves charge
    std::map<size_t,int> targetByZ;                            // species Z -> TARGET valence count N_val-q
    for (size_t i=0;i<atoms.size();i++)
        targetByZ[atoms[i].first] = atoms[i].second - q[i];
    return targetByZ;
}

// Shubnikov S3: the per-atom collinear labels, resolved by the SEED's own species rule (SeedCD's ctor):
// target = the ionic entry, else the neutral valence count; magnetic iff a spin pair exists AT that
// target (a stripped cation, target<=0, is non-magnetic).  Label = the atom's flip bit.
std::vector<int> MagneticDecoration(const Structure* st, const std::string& functional,
                                    const std::map<size_t,int>& ionicNvalByZ)
{
    assert(st);
    const std::string db = "atomic_valence_densities.json";
    std::map<size_t,bool> magneticByZ;
    std::vector<int> spins;
    st->ForEachSite([&](int iZ, const rvec3_t&, bool spinFlip)
    {
        const size_t Z = size_t(iZ);
        if (!magneticByZ.count(Z))
        {
            auto ti = ionicNvalByZ.find(Z);
            const int target = (ti!=ionicNvalByZ.end())
                             ? ti->second
                             : (int)std::lround(GetAtomicDensity(iZ, functional, db).Charge());
            magneticByZ[Z] = target>0 && HasAtomicSpinPair(iZ, functional, db, target);
        }
        spins.push_back(magneticByZ.at(Z) ? (spinFlip ? -1 : +1) : 0);
    });
    return spins;
}

} //namespace
