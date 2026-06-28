// File: ChargeDensity/Imp/Seed.C  SCF seed-density factory implementation (Phase 0: CoreGuess / Uniform).
module;
#include <cassert>
#include <type_traits>
#include <cstddef>

module qchem.ChargeDensity.Seed;
import qchem.ChargeDensity.Factory;   // IrrepCD_Factory<T>
import qchem.ChargeDensity.Types;     // tobs_t<T>
import qchem.Blaze;                    // blazem::zeroH, hmat_t
import qchem.Types;                    // dcmplx

namespace qchem::ChargeDensity
{

template <class T> tDM_CD<T>* MakeSeedDensity(SeedStrategy s, const BasisSet::BasisSet<T>* bs,
                                              const Structure* /*st -- Phases 1-3*/, const ElectronConfiguration* ec)
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
    case SeedStrategy::IonicSAD:
        assert(false && "SAD / IonicSAD seeds land in SCFSeedingPlan Phases 1-3");
        return nullptr;

    default:
        assert(false && "unknown SeedStrategy");
        return nullptr;
    }
}

template tDM_CD<double>* MakeSeedDensity<double>(SeedStrategy, const BasisSet::BasisSet<double>*,
                                                 const Structure*, const ElectronConfiguration*);
template tDM_CD<dcmplx>* MakeSeedDensity<dcmplx>(SeedStrategy, const BasisSet::BasisSet<dcmplx>*,
                                                 const Structure*, const ElectronConfiguration*);

} //namespace
