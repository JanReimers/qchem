// File: BasisSet/Gaussian/Lattice/Imp/RasterLadder.C  The GPW raster ladder (V-BS-RASTER).  The arithmetic is the evaluator's, step for
// step (its ctor's density-grid build, PairCap, MakeLevel, BuildLevels), so the two agree to the bit -- UTLattice_BS RasterLadder.*
// compares them on Si Γ, Si k222, MnO and an O2 box.
module;
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <vector>
module qchem.BasisSet.Gaussian.Lattice.RasterLadder;
import qchem.BasisSet.Orbital_1E_IBS;   // Real_OIBS (the molecular basis's one orbital block)

namespace qchem::BasisSet::Gaussian
{

namespace
{
//! The top completion rung's gate: a rung is appended when the reference grid sits below RelCutoffSafety() x kRungGateC x alpha_max
//! (the C=8-era measured calibration, deliberately decoupled from cutoffFactor -- see BuildLevels in Imp/GPW_Evaluator.C).
constexpr double kRungGateC=8.0;
}

GpwRasterLadder::GpwRasterLadder(const UnitCell& cell, const GaussianSharpness& basis, double densityEcut, double cutoffFactor,
                                 double ladderFactor, PlaneWave::RasterFields rasterFields, const GPWTolerances& tol)
    : itsInputs{densityEcut, cutoffFactor, ladderFactor, rasterFields, tol}
    , itsAlphaMax(basis.MaxExponent()), itsAlphaMin(basis.MinExponent()), itsRelCutoffSafety(basis.RelCutoffSafety())
    , itsLadderFactor(ladderFactor), itsTol(tol)
{
    if (rasterFields==PlaneWave::RasterFields::HartreeOnly) itsRelFieldSharp=itsTol.relFieldSharp*itsAlphaMax;
    if (densityEcut==0.0) return;                       // the DFT tier is off: no grid, no ladder
    const double aMax  = itsAlphaMax;
    const double floor = cutoffFactor * aMax;           // cutoffFactor*alpha_max (DENSITY scale)
    double ecut = densityEcut;
    if (densityEcut < 0.0) ecut = floor;                // AUTOMATIC: the basis-derived floor
    else if (densityEcut < floor)                       // EXPLICIT but under-resolved: warn, honour as given
        std::cerr << "[GPW] WARNING: densityEcut=" << densityEcut << " < " << cutoffFactor
                  << "*alpha_max=" << floor << " (basis alpha_max=" << aMax << "): the density grid "
                  << "under-resolves the basis -- charge will leak off-grid (integral rho_grid < N). Prefer "
                  << "densityEcut>=" << floor << ", or the automatic default (densityEcut<0)." << std::endl;
    // L0's raster: the density demand at the exponent L0 serves -- every pair when it is the finest level (an explicit mgridEcuts list,
    // or no completion rung above it), else the relative rule's cap -- raised for the uniform XC's own demand.
    const bool rung = itsTol.mgridEcuts.empty() && ecut < itsRelCutoffSafety*kRungGateC*aMax;
    itsGrid=MakeLevel(ReciprocalLattice(cell.MakeReciprocalCell()), ecut, PairCap(ecut, ecut, !rung), std::max(itsTol.rasterX, itsTol.xcRasterX));
    itsLevels=MakeLevels(itsGrid);
}

const RasterLevels& GpwRasterLadder::Levels() const
{
    if (!itsGrid) throw std::logic_error("GpwRasterLadder::Levels: the DFT tier is off (densityEcut 0) -- there is no ladder");
    return itsLevels;
}

double GpwRasterLadder::PairCap(double ecut, double ecutRef, bool finest) const
{
    const double amax=itsAlphaMax;
    const double beta=itsRelFieldSharp>=0.0 ? itsRelFieldSharp : itsTol.fieldSharp*amax;   // the relative rule's field floor
    const double pmax=std::max(2.0*amax, beta);                                              // the tightest p_eff of any pair
    if (finest) return pmax;
    const double cap = itsTol.relCutoff>0.0 ? ecut/itsTol.relCutoff
                                            : ecut*2.0*amax/(itsRelCutoffSafety*ecutRef);
    return std::min(pmax, cap);
}

std::shared_ptr<const PlaneWave::PW_Grid_Evaluator>
GpwRasterLadder::MakeLevel(const ReciprocalLattice& recip, double ecut, double pCap, double x) const
{
    assert(x>=0.0 && pCap>0.0);
    return std::make_shared<const PlaneWave::PW_Grid_Evaluator>(recip, rvec3_t(0,0,0), ecut, PlaneWave::NyquistRaster{std::sqrt(x*pCap)});
}

RasterLevels GpwRasterLadder::MakeLevels(std::shared_ptr<const PlaneWave::PW_Grid_Evaluator> grid) const
{
    if (!grid) throw std::invalid_argument("GpwRasterLadder::MakeLevels: no fine grid");
    RasterLevels L;
    if (!itsTol.mgridEcuts.empty())                     // GRID-MATCHING OVERRIDE: the explicit list replaces the progression + the rung
    {
        L.levels.push_back(grid);
        for (const double e : itsTol.mgridEcuts)
        {
            assert(e>0.0 && "GPWTolerances::mgridEcuts: sub-level cutoffs must be positive");
            if (e>=grid->Ecut())
            {
                std::cerr<<"[GPW] mgridEcuts: skipping sub-level "<<e<<" >= reference "<<grid->Ecut()<<std::endl;
                continue;
            }
            L.levels.push_back(MakeLevel(grid->Recip(), e, PairCap(e, grid->Ecut(), false), itsTol.rasterX));
        }
        L.nBase=L.levels.size();
    }
    else
    {
        const double efine=grid->Ecut();
        const double amax=itsAlphaMax, amin=itsAlphaMin;
        const double ecoarse=efine*amin/amax;           // resolves the most-diffuse pair product (exponent 2*amin)
        L.levels.push_back(grid);
        const double f=itsLadderFactor>1.0 ? itsLadderFactor : 4.0;
        double e=efine;
        while (e/f>=ecoarse)                            // factor-f coarsening down to the diffuse floor
        {
            e/=f;
            auto g=MakeLevel(grid->Recip(), e, PairCap(e, efine, false), itsTol.rasterX);
            static constexpr int kMinLevelN=3;          // the DEGENERACY FLOOR: a cell too small to hold a coarser grid
            const ivec3_t Ng=g->FFTGrid();
            if (std::min(std::min(Ng.x,Ng.y),Ng.z) < kMinLevelN) break;
            L.levels.push_back(g);
        }
        L.nBase=L.levels.size();
        if (efine < itsRelCutoffSafety*kRungGateC*amax) // the TOP COMPLETION RUNG
        {
            const double e2=itsRelCutoffSafety*efine;
            L.levels.push_back(MakeLevel(grid->Recip(), e2, PairCap(e2, efine, true), itsTol.rasterX));
        }
    }
    for (const auto& g : L.levels) { L.N.push_back(g->FFTGrid()); L.ecut.push_back(g->Ecut()); }
    return L;
}

RasterLevels GpwRasterLadder::MakeVlocLongLevels(double beta, double lnE) const
{
    const RasterLevels& D=Levels();
    RasterLevels L=D;
    // The custom top level: the sharpest pair's harmonic requirement (2 alpha_max against beta) -- appended LAST so ecut[0] stays the
    // resolution reference.  Skipped when an existing level (the completion rung) already satisfies it.  The requirement SATURATES at
    // 2 lnE beta for a hard basis (the field's own bandwidth is the ceiling), so the ball never runs away.
    const double p2=2.0*itsAlphaMax;
    const double ecutTop=2.0*lnE*p2*beta/(p2+beta);
    double emax=0.0;
    for (double e : L.ecut) emax=std::max(emax,e);
    const double xV=std::max(itsTol.vlocRasterX, itsTol.rasterX);
    if (ecutTop>emax)
    {
        auto g=MakeLevel(itsGrid->Recip(), ecutTop, PairCap(ecutTop, L.ecut[0], true), xV);
        L.levels.push_back(g); L.N.push_back(g->FFTGrid()); L.ecut.push_back(g->Ecut());
    }
    // V_long's OWN rasters (C2 step 3): when its demand exceeds the density's, every level of this private ladder is re-sampled at it
    // (the same balls, so the harmonic routing is unchanged; a static, once-per-geometry integrate-back).
    if (itsTol.vlocRasterX>itsTol.rasterX)
    {
        size_t Lf=0;
        for (size_t l=1;l<L.ecut.size();l++) if (L.ecut[l]>L.ecut[Lf]) Lf=l;
        const size_t nOwn = ecutTop>emax ? L.levels.size()-1 : L.levels.size();   // the top level was built at xV already
        for (size_t l=0;l<nOwn;l++)
        {
            L.levels[l]=MakeLevel(itsGrid->Recip(), L.ecut[l], PairCap(L.ecut[l], L.ecut[0], l==Lf), xV);
            L.N[l]=L.levels[l]->FFTGrid();
        }
    }
    return L;
}

std::shared_ptr<const GpwRasterLadder> MakeRasterLadder(const UnitCell& cell, const BasisSet::Real_BS& mol, const RasterLadderInputs& in)
{
    const BasisSet::Real_OIBS* only=nullptr;     // the single orbital block of the raw molecular basis (the evaluator's itsOrb)
    for (auto ibs : mol.Iterate<BasisSet::Real_OIBS>())
    {
        if (only) throw std::logic_error("MakeRasterLadder: a molecular basis with more than one orbital block");
        only=ibs;
    }
    const auto* sharp=dynamic_cast<const GaussianSharpness*>(only);   // abstract -> abstract
    if (!sharp) throw std::invalid_argument("MakeRasterLadder: the molecular basis has no Gaussian orbital block to size a raster from");
    return std::make_shared<const GpwRasterLadder>(cell, *sharp, in.densityEcut, in.cutoffFactor, in.ladderFactor, in.rasterFields, in.tol);
}

} // namespace qchem::BasisSet::Gaussian
