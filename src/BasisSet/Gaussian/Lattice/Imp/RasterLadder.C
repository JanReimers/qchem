// File: BasisSet/Gaussian/Lattice/Imp/RasterLadder.C  The GPW raster ladder (V-BS-RASTER).  The arithmetic is the evaluator's, step for
// step (its ctor's density-grid build, PairCap, MakeLevel, BuildLevels), so the two agree to the bit -- UTLattice_BS RasterLadder.*
// compares them on Si Γ, Si k222, MnO and an O2 box.
module;
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <memory>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
module qchem.BasisSet.Gaussian.Lattice.RasterLadder;
import qchem.BasisSet.Orbital_1E_IBS;   // Real_OIBS (the molecular basis's one orbital block)
import qchem.Reporting;                 // the grids facts and audit rows, written cursor-relative
import qchem.Symmetry.Lattice_3D.Fold;  // FoldGVectors (the grids report's {G}-star column)
import qchem.Matrix3D;                  // Transpose (the reciprocal face of a direct op)
import qchem.Vector3D;                  // norm (the raster Nyquist radius, |G|)
import qchem.Math;                      // Pi

namespace qchem::BasisSet::Gaussian
{

namespace
{
//! The top completion rung's gate: a rung is appended when the reference grid sits below RelCutoffSafety() x kRungGateC x alpha_max
//! (the C=8-era measured calibration, deliberately decoupled from cutoffFactor -- see BuildLevels in Imp/GPW_Evaluator.C).
constexpr double kRungGateC=8.0;

// THE GRID AUDIT (user 2026-10-25: "the max effective exponent and the ideal max G for every relevant H term ... quickly see
// that all grids are correct").  A term that integrates (or collocates) Gaussian PAIRS on a raster has TWO requirements:
//   BALL   -- its G-space sum must reach where the integrand's transform dies: a pair p against a field of Gaussian exponent
//             beta sums e^{-G^2 (1/4p + 1/4beta)}, the harmonic exponent p_h = p beta/(p+beta) (beta = 0: no field, p_h = p);
//             the tail at the ball edge is e^{-G_b^2/4p_h}.
//   RASTER -- the pair is SAMPLED in real space, so its spectrum beyond the raster's Nyquist G_N aliases back: the first
//             image sits at 2 G_N, relative amplitude e^{-(2G_N)^2/4p} = e^{-G_N^2/p}.  The ball's own raster puts G_N at the ball
//             edge (G_N^2 = 2 E_cut), and a raster demand x raises it to G_N^2 >= x p.  G_N here = the raster's inscribed Nyquist radius min_i pi N_i/|a_i|.
// The tightest pair (p = 2 alpha_max) is routed to the finest level, so the finest raster is the one that must pass.
// C2 step 3: the raster is SIZED to its demand x (GPWTolerances::rasterX ..., derived from the accuracy target), so the verdict is the
// demand's -- x = G_N^2/p >= xDemand, an ASSERTION (a level built short of it is a construction defect and throws).  The alias
// column stays the measured fact (the oracle gates: G_N^2/p = 2 is 4% on a tight pair's density and 1.5 Ha on its V_loc-long
// diagonal; 8 is <= 4e-5 / 3e-4 Ha -- UTLattice_BS GPW.TightPrimitive*).
double RasterNyquist(const ivec3_t& N, const UnitCell& A)
{
    const double a1=norm(A.ToCartesian(rvec3_t(1,0,0))), a2=norm(A.ToCartesian(rvec3_t(0,1,0))), a3=norm(A.ToCartesian(rvec3_t(0,0,1)));
    return std::min({Pi*N.x/a1, Pi*N.y/a2, Pi*N.z/a3});
}
//! The grid-audit ROW of one raster: PURE (a function of the raster and its demand).  THROWS when the raster is short of its demand
//! (a level built short of GPWTolerances' raster demand).
qchem::report::json GridAudit(const std::string& term, double p, double beta, double ecutBall, const ivec3_t& N, const UnitCell& A,
                              double xDemand, bool judgeBall, double ballTol=1e-6)
{
    namespace rpt = qchem::report;
    const double ph   = beta>0.0 ? p*beta/(p+beta) : p;
    const double Gb   = std::sqrt(2.0*ecutBall);
    const double tail = std::exp(-Gb*Gb/(4.0*ph));
    const double GN   = RasterNyquist(N, A);
    const double alias= std::exp(-GN*GN/p);
    if (GN*GN/p < xDemand*(1.0-1e-9))
        throw std::logic_error("GPW grid audit: "+term+" -- the finest raster delivers x = G_N^2/p = "+std::to_string(GN*GN/p)
                               +" below its demand "+std::to_string(xDemand)+" (a level was built short of GPWTolerances' raster demand)");
    const bool   ok   = !judgeBall || tail<=ballTol*(1.0+1e-9);   // the raster meets its demand (asserted above); the ball by the tolerance it was BUILT to
    rpt::json r;
    r["term"]        = term;
    r["pPair"]       = p;                                   // the tightest pair sampled on this raster (2 alpha_max)
    r["betaField"]   = beta;                                // the field's own Gaussian exponent (0 = none)
    r["ballEcut"]    = ecutBall;
    r["ballGneed"]   = std::sqrt(4.0*ph*std::log(1.0/ballTol));
    r["ballTol"]     = ballTol;
    r["ballGhave"]   = Gb;
    r["ballTail"]    = tail;
    r["ballJudged"]  = judgeBall;
    r["rasterN"]     = { N.x, N.y, N.z };
    r["rasterGneed"] = std::sqrt(p*xDemand);
    r["rasterGhave"] = GN;
    r["xDemand"]     = xDemand;
    r["x"]           = GN*GN/p;
    r["alias"]       = alias;
    r["status"]      = ok ? "OK" : "UNDER";
    return r;
}
//! File the row at <cursor>.<key>, and say UNDER loudly when it fails.  The ladder reports ONCE per point (it is shared), so there is
//! nothing to de-duplicate (the per-k repeat the old EmitAt dedup absorbed is gone with the per-block ladders).
void FileGridAudit(const std::string& key, const qchem::report::json& r)
{
    namespace rpt = qchem::report;
    rpt::Set(key, r);
    if (r["status"]=="OK") return;
    const auto N=r["rasterN"];
    std::ostringstream os;
    os<<"[grid audit] UNDER: "<<r["term"].get<std::string>()<<" -- tightest pair p="<<r["pPair"].get<double>()<<" on raster "
      <<N[0].get<int>()<<"x"<<N[1].get<int>()<<"x"<<N[2].get<int>()<<" (G_N="<<r["rasterGhave"].get<double>()
      <<", alias e^{-G_N^2/p}="<<r["alias"].get<double>()<<"): ball tail "<<r["ballTail"].get<double>()
      <<" (ball tolerance "<<r["ballTol"].get<double>()<<")";
    std::cout<<os.str()<<std::endl;
    rpt::Log(os.str());
}
size_t FinestLevel(const RasterLevels& L)
{
    size_t Lf=0;
    for (size_t l=1;l<L.ecut.size();l++) if (L.ecut[l]>L.ecut[Lf]) Lf=l;
    return Lf;
}
//! The reciprocal {U|τ} face of the imposed direct ops (U = Wᵀ, the G-index scatter map) -- what the {G}-star column folds under.
std::vector<Symmetry::Lattice_3D::SymOp> StarOps(const std::vector<Symmetry::Lattice_3D::DirectOp>& ops)
{
    std::vector<Symmetry::Lattice_3D::SymOp> sops;
    for (const auto& op : ops) sops.push_back({Transpose(op.W), op.tau});
    return sops;
}
}

GpwRasterLadder::GpwRasterLadder(const UnitCell& cell, const GaussianSharpness& basis, double densityEcut, double cutoffFactor,
                                 double ladderFactor, PlaneWave::RasterFields rasterFields, const GPWTolerances& tol, double ppSharpness)
    : itsInputs{densityEcut, cutoffFactor, ladderFactor, rasterFields, tol}
    , itsCell(cell), itsPPSharpness(ppSharpness)
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
    itsGrid=MakeLevel(ReciprocalLattice(cell.MakeReciprocal()), ecut, PairCap(ecut, ecut, !rung), std::max(itsTol.rasterX, itsTol.xcRasterX));
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

std::shared_ptr<const GpwRasterLadder> MakeRasterLadder(const UnitCell& cell, const BasisSet::Real_BS& mol, const RasterLadderInputs& in,
                                                        double ppSharpness)
{
    const BasisSet::Real_OIBS* only=nullptr;     // the single orbital block of the raw molecular basis (the evaluator's itsOrb)
    for (auto ibs : mol.Iterate<BasisSet::Real_OIBS>())
    {
        if (only) throw std::logic_error("MakeRasterLadder: a molecular basis with more than one orbital block");
        only=ibs;
    }
    const auto* sharp=dynamic_cast<const GaussianSharpness*>(only);   // abstract -> abstract
    if (!sharp) throw std::invalid_argument("MakeRasterLadder: the molecular basis has no Gaussian orbital block to size a raster from");
    return std::make_shared<const GpwRasterLadder>(cell, *sharp, in.densityEcut, in.cutoffFactor, in.ladderFactor, in.rasterFields, in.tol,
                                                   ppSharpness);
}

double LocalPPSharpness(const Structure& cl, const SpeciesRadialField& loc)
{
    // V_loc decays in G as e^{-G^2 rloc^2/2}: its effective Gaussian exponent is beta = 1/(2 rloc^2), the short-range Gaussian's alpha.
    double beta=0.0;
    if (const auto* gauss=dynamic_cast<const SpeciesRadialField_Gaussian*>(&loc))   // abstract -> abstract
        for (const Site& a : cl.Sites())
        {
            const auto terms=gauss->AsGaussians(a.Z,FieldRange::Short);
            if (!terms.empty()) beta=std::max(beta, terms[0].alpha);   // alpha = 1/(2 rloc^2)
        }
    return beta;
}

namespace
{
double VlocEps(const GPWTolerances& t) {return t.vlocEps>0.0 ? t.vlocEps : kVlocEpsNoTarget;}   // MakeLocalPPLong's rule
//! V_loc long's audit row: its finest private level at the tightest pair (p = 2 alpha_max) against the PP field beta, the ball judged
//! by the eps it was routed to.
qchem::report::json VlocLongRow(const RasterLevels& V, double p2, double beta, const UnitCell& cell, const GPWTolerances& t)
{
    const size_t Lf=FinestLevel(V);
    return GridAudit("V_loc long (integrate-back)", p2, beta, V.ecut[Lf], V.N[Lf], cell, std::max(t.vlocRasterX, t.rasterX),
                     /*judgeBall*/true, VlocEps(t));
}
//! The DENSITY raster's audit row: the finest density level at p = 2 alpha_max, against the routing field.  Only the RASTER is judged
//! here; the BALL is not "suppressed by 1/G^2" -- it enters the energy through the G-space quadratures (Hartree, V_loc long), and the
//! raster pilot judges it (SolidCalculation's CalibrateRasters: Eee at fixed D against twice the ball, doubling while over eps).  XC is
//! not where it enters: with the cusp-deficit feed (the default) XC reads rho[D] at the fixed point; linear D-mixing and GDM feed rho[D].
qchem::report::json DensityRow(const RasterLevels& L, double relFieldSharp, double aMax, const UnitCell& cell, double rasterX)
{
    const size_t Lf=FinestLevel(L);
    const double beta = relFieldSharp>=0.0 ? relFieldSharp : (2.0/3.0)*aMax;   // the routing field
    return GridAudit("density collocation (Hartree + XC raster)", 2.0*aMax, beta, L.ecut[Lf], L.N[Lf], cell, rasterX, /*judgeBall*/false);
}
}

std::vector<RasterVerdict> GpwRasterLadder::DeliveredRasterVerdicts() const
{
    std::vector<RasterVerdict> v;
    if (!itsGrid) return v;
    auto verdict=[](const char* audit, const qchem::report::json& r)
    { return RasterVerdict{audit, r["status"].get<std::string>(), r["alias"].get<double>()}; };
    v.push_back(verdict("auditDensity", DensityRow(itsLevels, itsRelFieldSharp, itsAlphaMax, itsCell, itsTol.rasterX)));
    if (itsPPSharpness>0.0)
        v.push_back(verdict("auditVlocLong", VlocLongRow(MakeVlocLongLevels(itsPPSharpness, -std::log(VlocEps(itsTol))), 2.0*itsAlphaMax,
                                                         itsPPSharpness, itsCell, itsTol)));
    return v;
}

void GpwRasterLadder::ReportGrids(const std::vector<Symmetry::Lattice_3D::DirectOp>& ops) const
{
    namespace rpt = qchem::report;
    rpt::Set("densityEcut",  itsInputs.cutoffFactor * itsAlphaMax);   // the auto density-grid floor Ecut
    rpt::Set("cutoffFactor", itsInputs.cutoffFactor);
    rpt::Set("raster",       { {"x", itsTol.rasterX}, {"vlocX", itsTol.vlocRasterX}, {"xcX", itsTol.xcRasterX} });   // the demands (0 = the ball's own)
    rpt::Set("routing",      itsRelFieldSharp>=0.0 ? "HartreeOnly" : "HartreeXC");   // field-sharpness routing policy
    rpt::Set("routingBeta",  itsRelFieldSharp>=0.0 ? itsRelFieldSharp : (2.0/3.0)*itsAlphaMax);
    if (!itsGrid) return;                                              // DFT tier off: no grids
    // The {G}-star column (see PrintGrids): star counts per level under the imposed ops, so the fold is VISIBLE in the grids table
    // (nG vs Gstars = the T1 reduction).  Empty ops = no column.
    const auto sops=StarOps(ops);
    const size_t Lf=FinestLevel(itsLevels);
    rpt::json ladder = rpt::json::array();
    for (size_t L = 0; L < itsLevels.levels.size(); ++L)
    {
        const PlaneWave::PW_Grid_Evaluator& lv = *itsLevels.levels[L];
        const ivec3_t N = lv.FFTGrid();
        // role: L==0 is the density/collocation reference (== FFT grid); the top rung completes the ladder above the local-PP
        // sub-ladder; everything between is a coarser multigrid level.
        const char* role = (L == 0) ? "reference" : (L >= itsLevels.nBase ? "rung" : "coarse");
        rpt::json row;
        row["level"] = (long)L;
        row["N"]     = { N.x, N.y, N.z };
        row["ecut"]  = lv.Ecut();
        row["nG"]    = (long)lv.size();
        {   // the alias exponent this level DELIVERS at the tightest pair the rule can route to it (C2 step 3: >= the demand)
            const double GN=RasterNyquist(N, itsCell);
            row["x"] = GN*GN/PairCap(lv.Ecut(), itsLevels.ecut[0], L==Lf);
        }
        if (!sops.empty())
            row["Gstars"] = (long)Symmetry::Lattice_3D::FoldGVectors(lv.Gs(), sops).Reps();
        row["role"]  = role;
        ladder.push_back(row);
    }
    rpt::Set("ladder", ladder);
    if (!sops.empty())
        rpt::Set("GstarFold", { { "ops", (long)sops.size() },
                                { "activeOn", "static local-PP sweeps (T1); per-iteration G fields unfolded" } });
    rpt::Set("localPP", { { "kappa", itsTol.localPPRelCutoff } });
    // the V_loc long routing tolerance that RUNS, and where it came from (MakeLocalPPLong's rule)
    rpt::json ve = { { "value",  VlocEps(itsTol) },
                     { "source", itsTol.VlocEpsDerived() ? "derived" : itsTol.vlocEps>0.0 ? "stated" : "no-target fallback" } };
    if (itsTol.VlocEpsDerived())   // derived = the rigorous per-pair map / the measured calibration (GridAccuracy::VlocEps)
    {
        ve["rigorous"]   = itsTol.vlocEpsRigorous;
        ve["etaRouting"] = itsTol.vlocRoutingEta;
    }
    rpt::Set("vlocEps", ve);
    FileGridAudit("auditDensity", DensityRow(itsLevels, itsRelFieldSharp, itsAlphaMax, itsCell, itsTol.rasterX));
    if (itsPPSharpness>0.0)            // V_long's row, a pure function of the ladder and the PP sharpness (no V_long run needed)
        FileGridAudit("auditVlocLong", VlocLongRow(MakeVlocLongLevels(itsPPSharpness, -std::log(VlocEps(itsTol))), 2.0*itsAlphaMax,
                                                   itsPPSharpness, itsCell, itsTol));
}

// GRID DIAGNOSTIC (doc/GPWPlan §0e): one line per stored grid, printed at run start so GPW's grids can be lined up against CP2K's
// &MGRID log (NGRIDS/CUTOFF/REL_CUTOFF + per-level N).  |G|min is the smallest NONZERO |G| (the cell scale 2*pi/a); |G|max the ball
// edge (sqrt(2*Ecut)).
std::ostream& PrintGrid(std::ostream& os, const std::string& tag, const PlaneWave::PW_Grid_Evaluator& g)
{
    double gmin=-1.0, gmax=0.0;
    for (const ivec3_t& m : g.Gs())
    {
        if (m.x==0 && m.y==0 && m.z==0) continue;
        const double gn=norm(g.GetGCartesian(m));
        if (gmin<0.0 || gn<gmin) gmin=gn;
        if (gn>gmax) gmax=gn;
    }
    const ivec3_t N=g.FFTGrid();
    return os<<"[GPW grid] "<<tag<<": N=("<<N.x<<","<<N.y<<","<<N.z<<") Ecut="<<g.Ecut()
             <<" Ha nG="<<g.size()<<" |G|min="<<gmin<<" |G|max="<<gmax<<"\n";
}

void GpwRasterLadder::PrintGrids(std::ostream& os, size_t nFunctions, const std::vector<Symmetry::Lattice_3D::DirectOp>& ops) const
{
    os<<"[GPW basis] n="<<nFunctions<<" alpha_min="<<itsAlphaMin<<" alpha_max="<<itsAlphaMax
      <<" cutoffFactor="<<itsInputs.cutoffFactor
      <<" raster x(density, V_long, XC L0)=("<<itsTol.rasterX<<", "<<itsTol.vlocRasterX<<", "<<itsTol.xcRasterX<<")"
      <<" routing="<<(itsRelFieldSharp>=0.0?"HartreeOnly beta=":"HartreeXC beta=2/3*alpha_max=")
      <<(itsRelFieldSharp>=0.0?itsRelFieldSharp:(2.0/3.0)*itsAlphaMax)<<" Ha-exp"
      <<" (auto density floor Ecut="<<itsInputs.cutoffFactor*itsAlphaMax<<" Ha)"<<std::endl;
    if (!itsGrid)
    {
        os<<"[GPW grid] DFT tier OFF (densityEcut=0): no grids"<<std::endl;
        return;
    }
    PrintGrid(os, "FFT rho<->G (density/collocation)", *itsGrid);
    // T1 {G}-STAR FOLD VISIBILITY (user 2026-08-15: "no hint anywhere of folding").  Under imposed ops each level's ball partitions
    // into stars; report the count + reduction beside nG.  HONESTY NOTE: today the fold is CONSUMED only by the static local-PP sweeps
    // (MakeLocalPP/MakeLocalPPLong T1 reduced evaluation); the per-iteration G-space consumers run unfolded -- doc/OpenWork.md item 5.
    const auto sops=StarOps(ops);
    for (size_t L=0;L<itsLevels.levels.size();L++)
    {
        std::string tag="ladder L="+std::to_string(L);
        if (!sops.empty())
        {
            const size_t stars=Symmetry::Lattice_3D::FoldGVectors(itsLevels.levels[L]->Gs(), sops).Reps();
            tag+=" stars="+std::to_string(stars);
        }
        if (L==0)                     tag+=" (== FFT grid; resolution reference)";
        else if (L>=itsLevels.nBase) tag+=" (top completion rung; local-PP uses base sub-ladder only)";
        PrintGrid(os, tag, *itsLevels.levels[L]);
    }
    os<<"[GPW grid] local-PP integration: FULL ladder L=0.."<<itsLevels.levels.size()-1
      <<" absolute REL_CUTOFF kappa="<<itsTol.localPPRelCutoff<<" Ha (e^{-kappa/2} pair tails)"<<std::endl;
    if (!sops.empty())
        os<<"[GPW grid] T1 {G}-star fold: "<<sops.size()<<" imposed ops; ACTIVE on the static local-PP"
            " sweeps (form factor at star reps); per-iteration G fields UNFOLDED (OpenWork item 5)"<<std::endl;
    os.flush();
}

RasterExponents GpwRasterLadder::DeliveredRasterExponents() const
{
    RasterExponents x;
    if (!itsGrid) return x;
    const size_t Lf=FinestLevel(itsLevels);
    const double Gf=RasterNyquist(itsLevels.N[Lf], itsCell), G0=RasterNyquist(itsLevels.N[0], itsCell);
    x.density =Gf*Gf/(2.0*itsAlphaMax);
    x.xcL0    =G0*G0/PairCap(itsLevels.ecut[0], itsLevels.ecut[0], Lf==0);
    if (itsPPSharpness>0.0)
    {
        const RasterLevels V=MakeVlocLongLevels(itsPPSharpness, -std::log(VlocEps(itsTol)));
        const double GN=RasterNyquist(V.N[FinestLevel(V)], itsCell);
        x.vlocLong=GN*GN/(2.0*itsAlphaMax);
    }
    return x;
}

} // namespace qchem::BasisSet::Gaussian
