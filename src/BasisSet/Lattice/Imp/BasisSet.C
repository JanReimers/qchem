// File: BasisSet/Lattice/Imp/BasisSet.C  Plane-wave basis-set container + factory implementation.
module;
#include <cassert>
#include <functional> // VetStageTrim's basis maker
#include <stdexcept>  // VetStageTrim's refusals
#include <cmath>     // lround (fractional k-point -> integer BZ-grid index)
#include <complex>   // std::norm (ShellConditioning's participation)
#include <limits>    // ShellConditioning: an empty block's lambda_min
#include <iostream>  // std::cout (the run-start GPW grid diagnostic)
#include <memory>    // std::shared_ptr / std::move (the GPW basis owns the molecular Gaussian basis)
#include <sstream>   // irrep label for the pre-flight basis.perIrrep rows
#include <optional>
#include <variant>   // the mixed child-slot visits (Step 3c-2)
#include <algorithm> // std::min (min singular value)
#include <cstdlib>   // std::getenv/std::atoi
#include <vector>    // std::vector (the k-block list + the space-group atom basis)
module qchem.BasisSet.Lattice.BasisSet;
import qchem.RunPolicy;   // theRunPolicy().StreamFold() -- the T3.2 fold, declared with the deviations (N5)
import qchem.BasisSet.Internal.BasisSetImp;   // BasisSetImp<dcmplx> (the generic list-of-IBS container)
import qchem.BasisSet.Gaussian.Lattice.GPW_IBS;     // Gaussian::GPW_IBS (the periodic-Gaussian block GPW_BasisSet owns)
import qchem.BasisSet.Orbital_1E_IBS;         // Real_OIBS (the molecular orbital block -- the stream-fold cross-cast)
import qchem.BasisSet.Gaussian.Lattice.LatticeSum1E;  // SetStreamSymmetryOps (the T3 route (b) stream fold, §6b)
import qchem.Reporting;                        // route the grid diagnostic into the run report when one is open
import qchem.Symmetry.Factory;                // BlochFactory (the Bloch irrep per k)
import qchem.Symmetry;                         // Spin, Irrep (the Bloch block identity for the pre-flight)
import qchem.Symmetry.Lattice_3D.SpaceGroup;   // SpaceGroup::Detect + AtomSite (crystal ops from the cell)
import qchem.Symmetry.Lattice_3D.BZReduction;  // ReduceToIBZ / IBZMesh (fold the MP mesh to the irreducible wedge)
import qchem.UnitCell;                          // UnitCell::GetCellMatrix / ToFractional (the space-group input)
import qchem.Structure;                         // Atom (itsZ / itsR -- the fractional basis for SpaceGroup)
import qchem.LASolver;                         // PivotedCholeskyDrops (the conditioning detector)
import qchem.BasisSet.AoShellSource;           // the shells a Bloch block is built from (VetStageTrim)
import qchem.Blaze;                            // blazem::eigen (grid-free overlap conditioning)
import qchem.Matrix3D;                          // Matrix3D (space-group cell matrix)
import qchem.Types;

namespace qchem::BasisSet::Lattice
{

// ONE plane-wave block per Brillouin-zone k-point: the basis ctor is the single place that enumerates k, so
// the framework's per-irrep loop (each k IS a Bloch irrep) becomes the BZ sum Sum_k w_k.  The KMesh carries
// the points + weights (uniform 1/Nk for an unreduced grid; symmetry-reduced points/weights will plug in
// here later).  N=(1,1,1) -> a single Gamma block.
PW_BasisSet::PW_BasisSet(const ::qchem::Lattice_3D& lat, double Ecut, const std::vector<int>& siteSpins)
{
    ReciprocalLattice recip=lat.Reciprocal();
    const ivec3_t N=lat.GetLimits();
    for (const auto& kp : lat.MakeKMesh())
    {
        ivec3_t ik(std::lround(kp.k.x*N.x), std::lround(kp.k.y*N.y), std::lround(kp.k.z*N.z));
        auto* pw=new PlaneWave::PlaneWave_IBS(recip, Symmetry::BlochFactory(N, ik, kp.weight), Ecut);
        Insert(pw);                                 // BasisSetImp takes ownership (no PP: the Vpseudo
                                                    // Hamiltonian term owns the pseudopotential model)
    }
    const bool magnetic=std::any_of(siteSpins.begin(), siteSpins.end(), [](int s){return s!=0;});
    if (!magnetic) return;
    if (GetNumIBS()!=1 || N!=ivec3_t(1,1,1))
        throw std::invalid_argument("PW_BasisSet: a magnetic imposition is Γ-only (one block, no Bloch phase per op)");
    // The group is resolved by the ONE imposition builder (closure-gated there), not here (V-BS-SYMMETRY S1).  The plane-wave
    // family imposes ONLY the spin pair so far: its density faces are cleared, so the composite density is not star-averaged
    // (exactly the behaviour before S3; the pair projection already makes each channel symmetric under the σ=None subgroup).
    // Its XC quadrature is the cell's plain mesh too (no Shubnikov fold: the δ route's (ρ,m) star-average is GPW's), so the
    // magnetic face is cleared with them and the imposition builds the plain mesh (V-BS-SYMMETRY S2).
    auto imp=std::make_shared<CrystalImposition>(ResolveImposition(lat, {.imposeSymmetry=true, .siteSpins=siteSpins}));
    imp->densityGroup=nullptr;
    imp->directDensity.clear();
    imp->magneticDirect.clear();
    SetImposition(std::move(imp));
    std::cout << "[symmetry] the plane-wave block imposes it by the (D_up, D_dn) pair projection" << std::endl;
}

Complex_BS* Factory(Type type, const ::qchem::Lattice_3D& lat, double Ecut, const std::vector<int>& siteSpins)
{
    assert(type==Type::PW);
    return new PW_BasisSet(lat, Ecut, siteSpins);
}

// GPW: one Bloch block of periodic Gaussians per Brillouin-zone k-point, built from the molecular basis over
// the cell's atoms -- the exact mirror of PW_BasisSet above (the ctor is the sole place that enumerates k, so
// the framework's per-irrep loop becomes the BZ sum Sum_k w_k).  N=(1,1,1) -> a single Gamma block.  The
// molecular basis is shared (shared_ptr) across the k-blocks.  Lattice images are eps-converged series
// summed inside the molecular seam (no radius exists); homeCellOnly is the finite "molecule in a box" MODE
// (image-free by definition -- every k-block identical), used by the finite==lattice gates.
// The k-mesh as (integer grid index, BZ weight) pairs -- either the FULL Monkhorst-Pack grid or, when
// p.imposeSymmetry, its irreducible wedge folded under the crystal point group (τ=0 symmorphic ops only -- the safe
// W-only guard; doc/GPWPlan1.md item 3).  IBZ needs the density symmetrized over each star to be exact; the
// caller is responsible for that (this only reduces which blocks are BUILT + carries the star weights).
struct KBlock { ivec3_t ik; double weight; };
// The run's imposed crystal ops (fold + density + Shubnikov + mesh invariance) are a qchem::CrystalImposition, resolved from the
// lattice and the options by ResolveImposition (qchem.Lattice_3D.Imposition; V-BS-SYMMETRY S1) -- not decided here.
//! The ladder settings a GPWParams states (what a handed-in ladder must have been built from).
static Gaussian::RasterLadderInputs LadderInputs(const GPWParams& p)
{
    return {.densityEcut=p.densityEcut, .cutoffFactor=p.cutoffFactor, .ladderFactor=p.ladderFactor, .rasterFields=p.rasterFields, .tol=p.tol};
}

static std::vector<KBlock> BuildKBlocks(const ::qchem::Lattice_3D& lat, const GPWParams& p,
                                        const CrystalImposition& cops)
{
    namespace SL = qchem::Symmetry::Lattice_3D;
    const std::vector<Matrix3D<double>>& ops = cops.recipFold;
    const rvec3_t kShift=p.kShift;
    const ivec3_t N=lat.GetLimits();
    std::vector<KBlock> blocks;
    if (!cops.policy.Imposes())
    {
        for (const auto& kp : lat.MakeKMesh(kShift))
            blocks.push_back({ivec3_t(std::lround(kp.k.x*N.x-kShift.x), std::lround(kp.k.y*N.y-kShift.y),
                                      std::lround(kp.k.z*N.z-kShift.z)), kp.weight});
        return blocks;
    }
    // IBZ: fold the MP grid under the (pre-computed) τ=0 reciprocal ops.  Every k-resolved BAND quantity obeys
    // f(Uk)=f(k) so the star weight is exact for the eigenvalue/occupation sum; the DENSITY needs the star
    // symmetrization (the composite ctor-injects the SAME ops, so a folded run is symmetrized automatically).
    SL::IBZMesh ibz = SL::ReduceToIBZ(N, kShift, ops);
    std::cout << "[IBZ] " << ibz.FullSize() << " k-points -> " << ibz.points.size()
              << " irreducible (point group |ops|=" << ops.size() << "); Σw=" << ibz.WeightSum() << std::endl;
    for (const auto& q : ibz.points) blocks.push_back({q.index, q.weight});
    return blocks;
}

GPW_BasisSet::GPW_BasisSet(const ::qchem::Lattice_3D& lat, std::shared_ptr<const BasisSet::Real_BS> mol,
                           const GPWParams& p)
{
    const rvec3_t kShift=p.kShift;
    const ivec3_t N=lat.GetLimits();
    // The run's imposition: HANDED IN (the facade resolves one per geometry point) or resolved here from the same options.  A
    // handed-in one must have been resolved from THESE options -- two disagreeing statements of the policy are refused, never
    // silently reconciled.
    const ImpositionOptions want{.imposeSymmetry=p.imposeSymmetry, .siteSpins=p.siteSpins,
                                 .freeMeshInvariant=p.freeMeshInvariant, .kShift=p.kShift};
    if (p.imposition && !(p.imposition->options==want))
        throw std::invalid_argument("GPW_BasisSet: GPWParams::imposition was resolved from other options than this GPWParams states "
                                    "(imposeSymmetry / siteSpins / freeMeshInvariant / kShift)");
    const std::shared_ptr<const CrystalImposition> imp = p.imposition ? p.imposition
                                                       : std::make_shared<const CrystalImposition>(ResolveImposition(lat, want));
    const CrystalImposition& ops = *imp;
    SetImposition(imp);                             // the density star ops + the spin-pair group, for the density / WF (S3)
    const std::vector<KBlock> blocks = BuildKBlocks(lat, p, ops);
    // The point's raster ladder (V-BS-RASTER): HANDED IN (the facade builds one per geometry point) or built here by the same builder.
    // Every block shares it -- the ladder's inputs are k-independent.  A handed-in one built from other settings is refused.
    if (p.ladder && !(p.ladder->Inputs()==LadderInputs(p)))
        throw std::invalid_argument("GPW_BasisSet: GPWParams::ladder was built from other settings than this GPWParams states "
                                    "(densityEcut / cutoffFactor / ladderFactor / rasterFields / tol)");
    const std::shared_ptr<const Gaussian::GpwRasterLadder> ladder = p.ladder ? p.ladder
                                                                 : Gaussian::MakeRasterLadder(lat.GetUnitCell(), *mol, LadderInputs(p));
    // T3 route (b) STREAM FOLD (doc/SymmetryUpgradePlan.md §6b/T3.2): on an IMPOSED Γ-ONLY run, fold the
    // shared molecular evaluator's (pair, R) collocation streams under the crystal ops -- reduced build +
    // replay (demand and per-iteration scatter/gather ÷ ~orbit factor); the existing per-iteration
    // SymmetrizeGMap/SymmetrizeRaster sites complete the group-average, so no raster commensurability is
    // needed.  Γ-only because a general-k block needs the little-group restriction + edge phases (T3.4);
    // multi-k IBZ runs keep full streams.
    // ARMED BY DEFAULT since 2026-08-19 (T3.5; OpenWork Step 2).  It was opt-in from 2026-08-03 because the
    // fold then imposed STRICTLY MORE than imposeSymmetry itself: the ρ star-average tolerates a
    // symmetry-BROKEN iterate D (it projects it pointwise), while the reduced replay read only
    // orbit-REPRESENTATIVE D elements -- i.e. it SAMPLED the orbit and thereby asserted D symmetric, which a
    // DEGENERATE OPEN SHELL at integer fill breaks permanently (measured: the imposed Si pseudo-atom-in-a-box
    // p² run flipped from the benign rotating-ρ mode into charge-transfer sloshing, ~0.26 Ha off).  The
    // replay now reads the ORBIT-PROJECTED D instead (NR_Evaluator::FoldProjectedD), and since collocation is
    // linear and equivariant that makes P ρ_red[D] == ρ[P D] == P ρ_full[D] for ANY iterate: the folded and
    // the unfolded imposed run compute the SAME density, so arming is a pure COST decision and needs no
    // auto-arm criterion.  GPW_STREAM_FOLD=0 is the opt-OUT, resolved through qchem::theRunPolicy() with the
    // other CP2K deviations (N5) rather than read here -- a gate that A/Bs it in one process calls
    // ReresolveRunPolicy() between the arms.
    if (ops.policy.Imposes() && blocks.size()==1
        && blocks[0].ik.x==0 && blocks[0].ik.y==0 && blocks[0].ik.z==0
        && kShift.x==0.0 && kShift.y==0.0 && kShift.z==0.0)
    {
        if (theRunPolicy().StreamFold())
        {
            // S3 guard: the stream fold asserts D ITSELF symmetric under its ops PER CHANNEL, and a
            // flip op relates D_up to D_dn -- so a magnetic imposition may fold streams only under the
            // σ=None (sublattice-preserving) subgroup.  Folding is merely reduced, never wrong.
            namespace SL = qchem::Symmetry::Lattice_3D;
            std::vector<SL::DirectOp> streamOps;
            if (ops.magneticDirect.empty()) streamOps = ops.directDensity;
            else for (const auto& op : ops.magneticDirect)
                if (op.sigma==SL::SpinAction::None) streamOps.push_back({op.W, op.tau});
            const BasisSet::Real_OIBS* orb=nullptr;
            for (auto ibs : const_cast<BasisSet::Real_BS&>(*mol).Iterate<BasisSet::Real_OIBS>()) { orb=ibs; break; }
            // THE NARROW FACE (ISP split 2026-09-08): arming the fold is all this does, so it names the
            // two-method StreamFoldable capability rather than the whole periodic aggregate.
            if (const auto* sf=dynamic_cast<const Gaussian::StreamFoldable*>(orb))
            {
                const size_t used=sf->SetStreamSymmetryOps(streamOps, lat.GetUnitCell());
                std::cout << "[stream fold] imposed Γ run: " << used << "/" << streamOps.size()
                          << " ops folded into the collocation streams" << std::endl;
            }
        }
    }
    for (const auto& kb : blocks)                                // the k-fold uses the linear reciprocal ops (+TR)
    {
        // Build the Bloch irrep WITH its BZ weight (star weight under IBZ) and the primary sym_t ctor -- the
        // weight carries the Sum_k w_k so the BZ-summed charge/energy are per-cell, not xNk.
        const sym_t irrep = Symmetry::BlochFactory(N, kb.ik, kb.weight, kShift);
        // THE FACTORY DECISION (doc/RealComplexPlan.md §1, Step 3c-3): a block is real ⇔ its irrep is
        // (TRIM, exact integer arithmetic) ∧ every Hamiltonian term preserves realness (the composition
        // root's fact, threaded in as GPWParams::hamPreservesReal).  One scalar-generic build serves
        // both alternatives; the typed Insert (Step 3b) files each under its own child slot.
        auto build=[&]<class U>()
        {
            auto* b=new Gaussian::tGPW_IBS<U>(lat.GetUnitCell(), irrep,
                                    mol, ladder, p.images,  // mol and the ladder shared across k-blocks
                                    ops.directDensity,    // {W|τ} = the IBZ raster star ops
                                    p.tol);               // the typed GPW tolerances (D-ENV step 5)
            Insert(b);
        };
        if (p.hamPreservesReal && irrep->IsReal()) build.template operator()<double>();
        else                                       build.template operator()<dcmplx>();
    }
    // GRID DIAGNOSTIC (doc/GPWPlan §0e): every stored grid, once per basis.  When a run report is open the orchestrator emits grids
    // EXPLICITLY (EmitGpwGrids) AFTER its conditioning pre-flight; only the legacy unbracketed console path stays here.
    if (!blocks.empty() && report::Depth() == 0) ladder->PrintGrids(std::cout, mol->GetNumFunctions(), ops.directDensity);
    // The Shubnikov PAIR projection: only where every op maps the one block onto itself with no Bloch phase -- a magnetic
    // imposition of a single Γ block.  That rule is the imposition's now (CrystalImposition::SpinPairOps), and the projection is
    // Symmetry::ProjectSpinPair over the block's own DenseRepresentation (V-BS-SYMMETRY S3).
}

// Compact irrep label for the report (symmetry symbol + spin arrow), via the Streamable Write().


std::shared_ptr<const Gaussian::GpwRasterLadder> MakeRasterLadder(const ::qchem::Lattice_3D& lat, std::shared_ptr<const BasisSet::Real_BS> mol,
                                                                  const GPWParams& p, double ppSharpness)
{
    if (!mol) throw std::invalid_argument("MakeRasterLadder: no molecular basis");
    return Gaussian::MakeRasterLadder(lat.GetUnitCell(), *mol, LadderInputs(p), ppSharpness);
}

void EmitGpwGrids(const Gaussian::GpwRasterLadder& ladder, const BasisSet::Real_BS& mol, const CrystalImposition& imp)
{
    if (report::Depth()==0) return;
    ladder.ReportGrids(imp.directDensity);
    // The lattice-sum economy readout (grids.latticeSums + the [lattice sums] console line): the numbers that JUMP when diffuse
    // functions are added -- the owner of the screens (the molecular basis's orbital block) reports them.
    for (auto ibs : mol.Iterate<BasisSet::Real_OIBS>())
        if (const auto* lat=dynamic_cast<const Gaussian::LatticeSum1E*>(ibs)) { lat->EmitLatticeSumReport(*imp.cell); break; }   // abstract -> abstract
}

// ------------------------------------------------------------------ the VET-STAGE TRIM (pin 22)

VetTrimResult VetStageTrim(const ::qchem::Lattice_3D& lat,
                           const std::function<std::shared_ptr<const BasisSet::Real_BS>(const Gaussian::ShellTrim&)>& make,
                           const GPWParams& p, double orthoTol)
{
    // The cell's sites, to name a shell's element from its centre.
    std::vector<std::pair<int,rvec3_t>> sites;
    lat.GetUnitCell().ForEachSite([&](int Z, const rvec3_t& R, bool){ sites.push_back({Z,R}); });
    auto elementAt=[&](const rvec3_t& c)
    {
        for (const auto& [Z,R] : sites)
        {
            const rvec3_t d=c-R;
            if (d.x*d.x+d.y*d.y+d.z*d.z < 1e-12) return Z;
        }
        throw std::logic_error("VetStageTrim: a basis shell whose centre is no site of the cell");
    };
    GPWParams full=p;
    full.imposeSymmetry=false;                 // the rank decision is the BASIS's: every k of the mesh, no fold
    full.siteSpins.clear();
    full.densityEcut=0.0;                      // the vet needs ANALYTIC overlaps only: no DFT tier, no grid ladder
    VetTrimResult r;
    // The smallest overlap eigenvalue over the k-mesh, and where -- the conditioning number a reader judges a
    // basis by (user 2026-09-28: the trim must show it BEFORE and AFTER).  One per pass; [0] = the untrimmed basis.
    struct LambdaMin { double value=1e300; std::string k; };
    std::vector<LambdaMin> lambdas;
    auto fmtLambda=[](const LambdaMin& l){ std::ostringstream os; os << l.value << " at k=" << l.k; return os.str(); };
    for (int pass=0; pass<64; ++pass)
    {
        r.mol=make(r.trim);
        GPW_BasisSet bs(lat, r.mol, full);
        const Complex_BS& cbs=bs;                  // the public face (GetNumIBS is protected on the Imp)
        // Per pass: every implicated (element, l), and the most diffuse shell of each (its exponents).
        struct Candidate { int Z; int l; rvec_t exponents; double alphaMin; };
        std::vector<Candidate> hits;
        size_t nDrops=0;
        LambdaMin lmin;
        for (size_t i=0;i<cbs.GetNumIBS();++i)
            std::visit([&](const auto& b)
            {
                const auto& S=b->Overlap();                                  // ANALYTIC, grid-free
                using U=typename std::decay_t<decltype(S)>::ElementType;
                {
                    rvec_t w; mat_t<U> Uw; blazem::eigen(S, w, Uw);          // ascending
                    if (w[0]<lmin.value) {lmin.value=w[0]; std::ostringstream os; os << b->GetSymmetry(); lmin.k=os.str();}
                }
                const std::vector<size_t> drops=qchem::PivotedCholeskyDrops<U>(S, orthoTol);
                if (drops.empty()) return;
                nDrops+=drops.size();
                const auto* src=dynamic_cast<const BasisSet::AoShellSource*>(b.get());   // abstract -> abstract
                if (!src) throw std::logic_error("VetStageTrim: a Bloch block that cannot name its shells (no AoShellSource)");
                const auto shells=src->GetAoShells();
                for (size_t idx : drops)
                    for (const auto& sh : shells)
                    {
                        if (idx<sh.offset || idx>=sh.offset+sh.nComponents()) continue;
                        const int l=sh.rep->L();
                        if (l>=2 && sh.nComponents()!=size_t(2*l+1))
                            throw std::runtime_error("VetStageTrim: a CARTESIAN l>=2 shell is implicated -- its l-2 "
                                "contaminants make 'the most diffuse shell of that l' the wrong cut (doc/Pins.md pin 22 (a)); "
                                "run the spherical view (GPW_SPHERICAL) for a vet-stage trim");
                        const int Z=elementAt(sh.center);
                        // The most diffuse shell of this (element, l) among ALL the block's shells.
                        Candidate best{Z,l,{},1e300};
                        for (const auto& t : shells)
                        {
                            if (t.rep->L()!=l || elementAt(t.center)!=Z || t.exponents.size()==0) continue;
                            double amin=t.exponents[0]; for (double a : t.exponents) amin=std::min(amin,a);
                            if (amin<best.alphaMin) best={Z,l,t.exponents,amin};
                        }
                        if (best.exponents.size()) hits.push_back(best);
                        break;
                    }
            }, bs.GetChild(i));
        lambdas.push_back(lmin);
        if (nDrops==0)
        {
            // THE RECORD (the caller opened the scope): the verdict a reader judges the basis by -- what went, and lambda_min of
            // S over the k-mesh before and after (user 2026-11-11, GUI item 27: "what was removed? ... am I good to go?").
            qchem::report::Set("orthoTol", orthoTol);
            qchem::report::Set("lambdaMinBefore", lambdas.front().value);
            qchem::report::Set("lambdaMinBeforeAt", lambdas.front().k);
            qchem::report::Set("lambdaMinAfter", lambdas.back().value);
            qchem::report::Set("lambdaMinAfterAt", lambdas.back().k);
            qchem::report::Set("nShellsRemoved", (long)r.trim.shells.size());
            std::cout << "[basis trim] vet stage: ";
            if (r.trim.empty()) std::cout << "nothing to trim";
            else { std::cout << "trimmed "; r.trim.Write(std::cout); }
            std::cout << " -- no k-block of the full mesh drops an AO at orthoTol=" << orthoTol
                      << ";  min eig S over the k-mesh: ";
            if (lambdas.size()==1) std::cout << fmtLambda(lambdas[0]);
            else std::cout << "BEFORE " << fmtLambda(lambdas.front()) << "  ->  AFTER " << fmtLambda(lambdas.back());
            std::cout << std::endl;
            return r;
        }
        if (hits.empty()) throw std::logic_error("VetStageTrim: drops that map to no nameable shell");
        // THE cut of this pass: the globally most diffuse implicated shell.
        const Candidate* cut=&hits[0];
        for (const auto& h : hits) if (h.alphaMin<cut->alphaMin) cut=&h;
        if (r.trim.Removes(cut->Z, cut->l, cut->exponents))
            throw std::logic_error("VetStageTrim: the trimmed shell is still in the basis (the reader ignored the trim)");
        r.trim.shells.push_back({cut->Z, cut->l, cut->exponents});
        {
            qchem::report::Row row("passes");                      // one per pass: why it cut, and the shell it cut
            std::ostringstream shell; Gaussian::ShellTrim one; one.shells.push_back(r.trim.shells.back()); one.Write(shell);
            qchem::report::Set("lambdaMin", lmin.value);
            qchem::report::Set("at", lmin.k);
            qchem::report::Set("aoDrops", (long)nDrops);
            qchem::report::Set("removes", shell.str());
            qchem::report::Set("Z", (long)cut->Z);
            qchem::report::Set("l", (long)cut->l);
            qchem::report::Set("alphaMin", cut->alphaMin);
        }
        std::cout << "[basis trim] vet stage pass " << pass << ": min eig S " << fmtLambda(lmin) << "; " << nDrops
                  << " AO drop(s) over the k-mesh at orthoTol=" << orthoTol << " -> removing ";
        Gaussian::ShellTrim one; one.shells.push_back(r.trim.shells.back()); one.Write(std::cout);
        std::cout << " from EVERY site of that element, every k" << std::endl;
    }
    throw std::runtime_error("VetStageTrim: 64 passes without a clean basis -- something other than diffuse shells is near-dependent");
}

//--------------------------------------------------------------------------------------------- THE SHELL WHAT-IF
ShellConditioning::ShellConditioning(const ::qchem::Lattice_3D& lat, std::shared_ptr<const BasisSet::Real_BS> untrimmed, const GPWParams& p)
{
    std::vector<std::pair<int,rvec3_t>> sites;
    lat.GetUnitCell().ForEachSite([&](int Z, const rvec3_t& R, bool){ sites.push_back({Z,R}); });
    auto elementAt=[&](const rvec3_t& c)
    {
        for (const auto& [Z,R] : sites)
        {
            const rvec3_t d=c-R;
            if (d.x*d.x+d.y*d.y+d.z*d.z < 1e-12) return Z;
        }
        throw std::logic_error("ShellConditioning: a basis shell whose centre is no site of the cell");
    };
    auto shellIndex=[&](int Z, int l, const rvec_t& e)       // one entry per (element, l, exponents): a basis fact, pin 22
    {
        for (size_t i=0;i<itsShells.size();++i)
        {
            const Shell& s=itsShells[i];
            if (s.Z!=Z || s.l!=l || s.exponents.size()!=e.size()) continue;
            bool same=true; for (size_t a=0;a<e.size();++a) same=same && s.exponents[a]==e[a];
            if (same) return i;
        }
        itsShells.push_back({Z,l,e,0});
        return itsShells.size()-1;
    };
    GPWParams full=p;
    full.imposeSymmetry=false;                 // every k of the mesh, no fold: the rank question is the BASIS's (as VetStageTrim)
    full.siteSpins.clear();
    full.densityEcut=0.0;                      // analytic overlaps only
    GPW_BasisSet bs(lat, untrimmed, full);
    const Complex_BS& cbs=bs;
    for (size_t i=0;i<cbs.GetNumIBS();++i)
        std::visit([&](const auto& b)
        {
            std::ostringstream os; os << b->GetSymmetry();
            itsK.push_back(os.str());
            itsS.push_back(b->Overlap());                                    // the block's own scalar, as the vet stage factorizes it
            const auto* src=dynamic_cast<const BasisSet::AoShellSource*>(b.get());   // abstract -> abstract
            if (!src) throw std::logic_error("ShellConditioning: a Bloch block that cannot name its shells (no AoShellSource)");
            std::vector<std::vector<size_t>> aos(itsShells.size());
            for (const auto& sh : src->GetAoShells())
            {
                const size_t k=shellIndex(elementAt(sh.center), sh.rep->L(), sh.exponents);
                if (aos.size()<itsShells.size()) aos.resize(itsShells.size());
                for (size_t c=0;c<sh.nComponents();++c) aos[k].push_back(sh.offset+c);
            }
            itsAOs.push_back(std::move(aos));
        }, bs.GetChild(i));
    for (auto& aos : itsAOs) aos.resize(itsShells.size());     // a shell first met in a later block: empty in the earlier ones
    if (!itsAOs.empty()) for (size_t k=0;k<itsShells.size();++k) itsShells[k].nAO=itsAOs[0][k].size();
}

std::vector<bool> ShellConditioning::KeptBy(const Gaussian::ShellTrim& trim) const
{
    std::vector<bool> kept(itsShells.size(), true);
    for (size_t k=0;k<itsShells.size();++k) kept[k]=!trim.Removes(itsShells[k].Z, itsShells[k].l, itsShells[k].exponents);
    return kept;
}

namespace
{
//! The AOs of the \a kept shells of one block, ascending (the trimmed basis keeps the untrimmed order).
std::vector<size_t> KeptAOs(const std::vector<std::vector<size_t>>& aos, const std::vector<bool>& kept)
{
    std::vector<size_t> idx;
    for (size_t k=0;k<aos.size();++k) if (kept[k]) idx.insert(idx.end(), aos[k].begin(), aos[k].end());
    std::sort(idx.begin(), idx.end());
    return idx;
}
//! S restricted to \a idx: exactly the overlap of the basis without the other shells.
template <class M> M Restricted(const M& S, const std::vector<size_t>& idx)   // M: hmat_t<T> (an alias T cannot be deduced through)
{
    M R(idx.size());
    for (size_t r=0;r<idx.size();++r) for (size_t c=r;c<idx.size();++c) R(r,c)=S(idx[r],idx[c]);
    return R;
}
}

ShellConditioning::Spectrum ShellConditioning::Measure(size_t block, const std::vector<bool>& kept, double floor) const
{
    const std::vector<size_t> idx=KeptAOs(itsAOs[block], kept);
    Spectrum out{std::numeric_limits<double>::infinity(), 0.0, 0, idx.size()};   // no function at this k: nothing near-dependent
    if (idx.empty()) return out;
    std::visit([&](const auto& S)
    {
        using T=typename std::decay_t<decltype(S)>::ElementType;
        const auto R=Restricted(S, idx);
        rvec_t w; blazem::eigen(R, w);                                       // ascending
        out.lambdaMin=w[0];
        out.cond=w[w.size()-1]/w[0];
        out.aoDrops=qchem::PivotedCholeskyDrops<T>(R, floor).size();         // the vet stage's own test
    }, itsS[block]);
    return out;
}

ShellConditioning::Answer ShellConditioning::Analyze(const std::vector<bool>& kept, double floor) const
{
    if (kept.size()!=itsShells.size())
        throw std::invalid_argument("ShellConditioning::Analyze: "+std::to_string(kept.size())+" flags for "+std::to_string(itsShells.size())+" shells");
    if (!(floor>0.0)) throw std::invalid_argument("ShellConditioning::Analyze: the floor must be > 0 (an orthoTol)");
    Answer a;
    size_t worst=0;
    a.lambdaMin=std::numeric_limits<double>::infinity();
    for (size_t b=0;b<itsS.size();++b)
    {
        const Spectrum sp=Measure(b, kept, floor);
        const Block blk{itsK[b], sp.nFunctions, sp.lambdaMin, sp.cond, sp.aoDrops};
        if (blk.lambdaMin<a.lambdaMin) { a.lambdaMin=blk.lambdaMin; a.worstK=blk.k; worst=b; }
        a.cond=std::max(a.cond, blk.cond);
        a.aoDrops+=blk.aoDrops;
        a.nFunctions=blk.nFunctions;
        a.blocks.push_back(blk);
    }
    // the near-null space at the worst k, and each shell's share of it
    a.participation.assign(itsShells.size(), 0.0);
    if (!itsS.empty())
    {
        const std::vector<size_t> idx=KeptAOs(itsAOs[worst], kept);
        if (!idx.empty())
            std::visit([&](const auto& S)
            {
                using T=typename std::decay_t<decltype(S)>::ElementType;
                rvec_t w; mat_t<T> V; blazem::eigen(Restricted(S, idx), w, V);   // ascending, columns = eigenvectors
                size_t nv=0; while (nv<w.size() && w[nv]<floor) ++nv;
                const size_t use=std::max<size_t>(nv,1);
                std::vector<size_t> pos(S.rows(), size_t(-1));
                for (size_t r=0;r<idx.size();++r) pos[idx[r]]=r;
                for (size_t k=0;k<itsShells.size();++k)
                    if (kept[k])
                        for (size_t ao : itsAOs[worst][k])
                            for (size_t v=0;v<use;++v) a.participation[k]+=std::norm(V(pos[ao],v))/double(use);
            }, itsS[worst]);
    }
    // what toggling each shell ALONE does, worst over k
    a.lambdaMinToggled.assign(itsShells.size(), std::numeric_limits<double>::infinity());
    a.condToggled.assign(itsShells.size(), 0.0);
    a.aoDropsToggled.assign(itsShells.size(), 0);
    for (size_t k=0;k<itsShells.size();++k)
    {
        std::vector<bool> t=kept; t[k]=!t[k];
        for (size_t b=0;b<itsS.size();++b)
        {
            const Spectrum sp=Measure(b, t, floor);
            a.lambdaMinToggled[k]=std::min(a.lambdaMinToggled[k], sp.lambdaMin);
            a.condToggled[k]=std::max(a.condToggled[k], sp.cond);
            a.aoDropsToggled[k]+=sp.aoDrops;
        }
    }
    // THE RECORD, one level in (the caller opened the scope): the verdict, then a row per k and a row per shell.
    qchem::report::Set("floor", floor);
    qchem::report::Set("healthy", a.Healthy());
    qchem::report::Set("lambdaMin", a.lambdaMin);
    qchem::report::Set("worstK", a.worstK);
    qchem::report::Set("cond", a.cond);
    qchem::report::Set("aoDrops", (long)a.aoDrops);
    qchem::report::Set("nFunctions", (long)a.nFunctions);
    for (const Block& blk : a.blocks)
    {
        qchem::report::Row row("blocks");
        qchem::report::Set("k", blk.k);
        qchem::report::Set("nFunctions", (long)blk.nFunctions);
        qchem::report::Set("lambdaMin", blk.lambdaMin);
        qchem::report::Set("cond", blk.cond);
        qchem::report::Set("aoDrops", (long)blk.aoDrops);
    }
    for (size_t k=0;k<itsShells.size();++k)
    {
        const Shell& sh=itsShells[k];
        qchem::report::Row row("shells");
        qchem::report::Set("Z", (long)sh.Z);
        qchem::report::Set("l", (long)sh.l);
        std::vector<double> e; for (size_t i=0;i<sh.exponents.size();++i) e.push_back(sh.exponents[i]);
        qchem::report::Set("exponents", e);
        qchem::report::Set("nAO", (long)sh.nAO);
        qchem::report::Set("kept", bool(kept[k]));
        qchem::report::Set("lambdaMinToggled", a.lambdaMinToggled[k]);
        qchem::report::Set("condToggled", a.condToggled[k]);
        qchem::report::Set("aoDropsToggled", (long)a.aoDropsToggled[k]);
        qchem::report::Set("participation", a.participation[k]);
    }
    return a;
}

Complex_BS* GPWFactory(const ::qchem::Lattice_3D& lat, std::shared_ptr<const BasisSet::Real_BS> mol,
                       const GPWParams& p)
{
    return new GPW_BasisSet(lat, std::move(mol), p);
}
Complex_BS* GPWFactory(const ::qchem::Lattice_3D& lat, std::shared_ptr<const BasisSet::Real_BS> mol,
                       double densityEcut, rvec3_t kShift, Gaussian::CellImages images, double cutoffFactor)
{
    return new GPW_BasisSet(lat, std::move(mol),
                            GPWParams{densityEcut, cutoffFactor, images, kShift});
}

} //namespace
