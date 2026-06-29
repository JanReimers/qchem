// File: ChargeDensity/Imp/Seed.C  SCF seed-density factory implementation (Phase 0: CoreGuess / Uniform).
module;
#include <cassert>
#include <type_traits>
#include <cstddef>
#include <memory>
#include <vector>
#include <utility>
#include <algorithm>
#include <map>
#include <cmath>

module qchem.ChargeDensity.Seed;
import qchem.PeriodicTable;                   // thePeriodicTable -> GetElectronegativity (the ionic heuristic)
import qchem.ChargeDensity.Factory;          // IrrepCD_Factory<T>
import qchem.ChargeDensity.Types;            // tobs_t<T>
import qchem.ChargeDensity.CompositeFittedCD;// CompositeFittedCD (the molecular SAD seed, double only)
import qchem.ChargeDensity.FourierSeedCD;    // FourierSeedCD (the plane-wave SAD seed, dcmplx only)
import qchem.ChargeDensity.AtomicDensity;    // GetAtomicDensity, RadialDensity, RecentredAtomicDensity
import qchem.BasisSet.Band_FT_IBS;           // Band_FT_IBS (the G-space block for the PW seed)
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

template <class T> tChargeDensity<T>* MakeSeedDensity(SeedStrategy s, const BasisSet::BasisSet<T>* bs,
                                                      const Structure* st, const ElectronConfiguration* ec)
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
        const tobs_t<T>* block = (*bs)[0];
        const Irrep      irr   = bs->GetIrreps(Spin::None)[0];
        const size_t     n     = block->GetNumFunctions();
        const int        N     = ec->GetN(irr);
        hmat_t<T> D0=blazem::zeroH<T>(n);
        for (size_t i=0;i<n;i++) D0(i,i)=T(double(N)/double(n));
        return IrrepCD_Factory<T>(D0, block, irr);
    }

    case SeedStrategy::SAD:
    {
        // Superposition of neutral atomic densities (the molecular SAD seed, double/AO path).  Drop each
        // atom's recentred radial density (LDA DB) into a CompositeFittedCD; FittedVee/FittedVxc consume it
        // (Coulomb projection + rho(r)).  The plane-wave (dcmplx, FT) face is Phase 2.
        if constexpr (std::is_same_v<T,double>)
        {
            assert(st && "SAD seed needs a Structure (atom Z + positions)");
            auto* cd = new CompositeFittedCD(st->GetNumElectrons());
            for (size_t i=0;i<st->GetNumAtoms();i++)
            {
                const Atom* a = (*st)[i];
                auto rad = std::make_shared<const RadialDensity>(GetAtomicDensity(a->itsZ));
                cd->Insert(std::make_shared<const RecentredAtomicDensity>(rad, a->itsR));
            }
            return cd;
        }
        else
        {
            // Plane-wave (FT) SAD: rho-tilde(G) = Sum_atoms rho_atom(|G|) e^{-iG.R}, assembled by the basis.
            assert(st && "SAD plane-wave seed needs a Structure");
            const auto* ftbs = dynamic_cast<const BasisSet::Band_FT_IBS*>((*bs)[0]);
            assert(ftbs && "SAD plane-wave seed needs a Band_FT_IBS (plane-wave) basis");
            return new FourierSeedCD(ftbs, st);
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
            const auto* ftbs = dynamic_cast<const BasisSet::Band_FT_IBS*>((*bs)[0]);
            assert(ftbs && "IonicSAD plane-wave seed needs a Band_FT_IBS (plane-wave) basis");

            std::vector<std::pair<int,int>> atoms;                      // {Z, N_val} per atom
            std::map<size_t,int> nvalByZ;
            for (size_t i=0;i<st->GetNumAtoms();i++)
            {
                int Z = (*st)[i]->itsZ;
                if (!nvalByZ.count(Z))
                    nvalByZ[Z] = (int)std::lround(GetAtomicDensity(Z,"LDA","atomic_valence_densities.json").Charge());
                atoms.emplace_back(Z, nvalByZ[Z]);
            }
            std::vector<int> q = IonicFormalCharges(atoms);             // Na+1, F-1; conserves charge
            std::map<size_t,double> scaleByZ;
            for (size_t i=0;i<atoms.size();i++)                         // species Z -> (N_val - q)/N_val
                scaleByZ[atoms[i].first] = double(atoms[i].second - q[i]) / double(atoms[i].second);
            return new FourierSeedCD(ftbs, st, "LDA", scaleByZ);
        }
        else
        {
            assert(false && "molecular (AO) IonicSAD is a later phase; use SAD or the plane-wave path");
            return nullptr;
        }
    }

    default:
        assert(false && "unknown SeedStrategy");
        return nullptr;
    }
}

template tChargeDensity<double>* MakeSeedDensity<double>(SeedStrategy, const BasisSet::BasisSet<double>*,
                                                        const Structure*, const ElectronConfiguration*);
template tChargeDensity<dcmplx>* MakeSeedDensity<dcmplx>(SeedStrategy, const BasisSet::BasisSet<dcmplx>*,
                                                        const Structure*, const ElectronConfiguration*);

} //namespace
