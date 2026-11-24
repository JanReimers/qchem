// File: Structure/Imp/MeshFactory.C  The structure mesh factory's bodies (V-MESH-FACTORY §8 step 2): the geometry builders MOVED VERBATIM
// from Imp/UnitCell.C (the periodic Becke grid, the uniform raster, the site-adapted imposed build) and Imp/Molecule.C (the atom-centred
// Becke grid), with the two UnitCell members rewritten as free MakeMesh overloads over `cell`.  Same points, same order, same weights.
module;
#include <cstdio>
#include <iomanip>
#include <cassert>
#include <cstdlib>
#include <stdexcept> // RequireSiteBlocks -- an atom-centred-mesh invariant that must survive NDEBUG
#include <iostream>
#include <ostream>
#include <string>
#include <vector>
#include <functional>
#include <exception> // std::exception_ptr (the threaded weight response)
#include <algorithm> // std::any_of (is the fold's op set magnetic?)
#include <sstream>   // the molecular frame's report string
#include <cmath>
#include <limits>    // PointGroupGeometryError: an op set outside D2h type is infinitely far
#include <utility>
#include <memory>    // the scoped-but-early-stopped Timed on the site-adapted angular build (1.1(b))
module qchem.Structure.MeshFactory;
import qchem.Diagnostics;                     // the ONE diagnostics registry (D-ENV)
import qchem.Parallel;       // WorkerThreads (QCHEM_OPENMP_THREADS: the Becke mesh build)
import qchem.Vector3D;       // norm(rvec3_t)
import qchem.Mesh.Product;   // ProductMesh, MakeRadial, MakeAngular (the single-centre template)
import qchem.Mesh.Builder;   // MeshBuilder (efficient incremental accumulation)
import qchem.Mesh.BeckePartition;   // BeckePartitionWeight / PeriodicBeckeWeight: the ONE partition (the build and its derivative)
import qchem.SymmetrizeMesh; // MakeInvariantAngularMesh + FoldPointsPeriodic + MakeInvariant (the §6a imposed builds)
import qchem.Symmetry.Group; // Symmetry::Frame, TrivialGroup, kAdmissionTolerance (what each mesh forces, pin 30)
import qchem.Reporting;      // the meshes' fields, one level into the caller's report scope (inert when no run is open)
import qchem.Symmetry.Lattice_3D.Fold;   // FlipFixedPointsPeriodic (a Shubnikov fold's m = 0 points)
import qchem.Matrix3D;                   // Invert (fractional mesh coordinates)

namespace qchem {

// ---- The periodic Becke integration mesh (mp.cellKind==UnitCellKind::Becke) -------------------------------
// The atom-centred fuzzy-Voronoi quadrature of MakeMolecularMesh (Becke 1988), extended to a periodic cell.
// Each cell atom a carries the single-centre radial x angular product grid; a point r on atom a's grid gets
// the Becke cell weight
//     w_a(r) = P_a0(r) / sum_{b,L} P_bL(r),    P_bL(r) = prod_{(c,M)!=(b,L)} s(mu),
// where the competitor set runs over ALL periodic images R_b + L of the cell atoms (L = lattice vectors).
// The lattice translates of the per-atom weighted grids form a partition of unity over the whole crystal
// (sum over ALL images of w == 1 pointwise), so for a cell-periodic integrand f
//     sum_{a in cell} integral w_a f d^3r  =  integral_cell f d^3r
// exactly -- the per-cell weighted grids ARE the cell quadrature.  Points are wrapped into the home cell
// on emission (the weight is lattice-invariant, so wrapping is free).
//
// The image series is an epsilon-CONVERGED series, never radius-truncated (GPWPlan durable pin: magnitude
// screening is the only truncation mechanism).  Competitors are gathered in expanding Chebyshev cell shells
// around the point's own cell until BOTH tails are provably below eps:
//   (i)  any unseen image's own cell function:  P  <=  s((D-dn)/(D+dn))            (opposite-side collinear
//        bound; dn = distance to the nearest seen competitor, D = shell distance floor), and
//   (ii) any unseen image's FACTOR on a retained P_b:  P_b-bound * (1 - s((d_b-D)/(d_b+D)))  <  eps
// -- both evaluated with the same s(mu) polynomial, so the screen is a magnitude bound with no radius
// anywhere.  Far radial-tail points whose own cell function is below eps are dropped up front (some image
// atom is far closer, so the partition assigns them ~0 weight): a diffuse tail costs almost nothing, which
// is the point of the Becke XC grid (doc/GPWPlan1.md).
namespace {


//! GPW_BECKE_COUNT=1: the per-point COST CENSUS of the partition (doc/OpenWork.md Step 3).  The
//! partition is the single largest CPU consumer in the code (~50% of cycles), and it is threaded by
//! default -- so it reads as a small WALL bucket and the ledger cannot say where inside it the time
//! goes.  Two call sites compete for the blame: the shell-convergence RETEST (O(|im|) per shell, and
//! |im| grows as s^3) and the final P-set DOUBLE LOOP (|P-set| x |im|).  This census separates them
//! by counting, which is exact and thread-independent, instead of by profiling a threaded loop.
//! Slot-indexed per point and reduced after the loop -- the same discipline as keep/kpt/kwt above --
//! so it needs no atomics and cannot perturb what it measures.
//! Only the INNER-loop counter is gated on \a census: it fires ~2.3e5 times per point (97% of the work),
//! where an ungated increment would be a real production cost.  The shell-retest and P-set counters fire
//! ~6.4e3 and ~1e2 times per point, measured at ~0.03% of the build, and stay ungated for legibility.
using qcMesh::BeckePointCost;   // the census (qcMesh::PeriodicBeckeWeight fills it)

//! \brief THE PERIODIC BECKE MESH's WEIGHT RESPONSE (qcMesh::MeshWeightResponse; lane B 2026-11-19, the XC term's grid response on the
//! Becke mesh).  Installed by MakePeriodicBeckeMesh, which alone knows how the weights were made: each kept point g of site a carries
//! w_g = w_q W_a(R_a + u_g), w_q its template weight, W_a the periodic partition (qcMesh::PeriodicBeckeWeight).  So for f held per point,
//! ∂(Σ_g f_g w_g)/∂X = Σ_g f_g w_q ∂W_a/∂X over every competitor image X (the point riding its owner: qcMesh::PeriodicBeckeWeightGradient),
//! folded onto the cell atoms for the gradient and contracted with X for the virial (the template offsets u held: they do not strain).
//! The partition is RE-FORMED at R_a + u_g, the same arithmetic the build ran, so the competitor and P-sets are the build's.
class PeriodicBeckeResponse : public qcMesh::MeshWeightResponse
{
public:
    PeriodicBeckeResponse(const CellGeometry& cell, std::vector<rvec3_t> R, double dPlane, int k, double eps)
        : itsCell(cell), itsR(std::move(R)), itsDPlane(dPlane), itsK(k), itsEps(eps) {}
    Result Pullback(const qcMesh::Mesh& m, const rvec_t& f) const override
    {
        const size_t natom=itsR.size();
        if (f.size()!=m.size()) throw std::invalid_argument("PeriodicBeckeResponse: the field is "+std::to_string(f.size())+" for "+std::to_string(m.size())+" points");
        if (m.OwnerOffsets().size()!=m.size() || m.NSites()!=natom)
            throw std::logic_error("PeriodicBeckeResponse: the mesh carries no owner offsets or not one site block per atom");
        qchem::report::Timed timed("force: becke weight response");
        const rvec3vec_t& U=m.OwnerOffsets();
        const rvec_t& W=m.Weights();
        Result out{rvec3vec_t(natom, rvec3_t(0,0,0)), rmat3d_t(0,0,0, 0,0,0, 0,0,0)};
        double V[3][3]={{0,0,0},{0,0,0},{0,0,0}};
        auto point=[&](size_t a, size_t g, rvec3vec_t& G, double (&Vt)[3][3])
        {
            if (f[g]==0.0) return;
            std::vector<qcMesh::BeckeImage> im; long target=-1; BeckePointCost pc;
            const rvec3_t r=itsR[a]+U[g];
            const double wP=qcMesh::PeriodicBeckeWeight(itsCell, itsR, itsDPlane, a, r, itsK, itsEps, im, target, pc, false);
            if (!(wP>0.0) || target<0) return;   // the build kept it, so this is unreachable short of a changed geometry
            // THE PLATEAUS (cost; measured 2026-11-19, the AD was the whole Becke force): where W is 1 or 0 to 1e-14 the partition's other
            // cell functions are below 1e-14 of the point's own (or its own below theirs), and W's derivative is of that order -- a point
            // deep in its own basin, the bulk of a radial grid, contributes nothing to its weight response.
            if (wP>1.0-1e-14 || wP<1e-14) return;
            const double s=f[g]*W[g]/wP;          // f_g w_q
            const std::vector<rvec3_t> dW=qcMesh::PeriodicBeckeWeightGradient(im, size_t(target), U[g], itsK);
            for (size_t i=0;i<im.size();i++)
            {
                const rvec3_t gi=s*dW[i];
                G[im[i].atom]+=gi;
                const double gg[3]={gi.x,gi.y,gi.z}, xx[3]={im[i].R.x,im[i].R.y,im[i].R.z};
                for (int x=0;x<3;x++) for (int y=0;y<3;y++) Vt[x][y]+=gg[x]*xx[y];
            }
        };
        std::vector<std::pair<size_t,size_t>> work;   // (site, point)
        for (size_t a=0;a<natom;a++) for (size_t g=m.SiteBegin(a); g<m.SiteEnd(a); g++) work.push_back({a,g});
#ifdef QCHEM_OPENMP
        const int nthreads=qchem::InnerThreads();
        if (nthreads>1)
        {
            std::exception_ptr firstEx;
            #pragma omp parallel num_threads(nthreads)
            {
                rvec3vec_t Gt(natom, rvec3_t(0,0,0));
                double Vt[3][3]={{0,0,0},{0,0,0},{0,0,0}};
                #pragma omp for schedule(dynamic, 16)
                for (size_t w=0; w<work.size(); w++)
                {
                    try { point(work[w].first, work[w].second, Gt, Vt); }
                    catch (...)
                    {
                        #pragma omp critical (becke_response_throw)
                        if (!firstEx) firstEx=std::current_exception();
                    }
                }
                #pragma omp critical (becke_response_reduce)
                {
                    for (size_t a=0;a<natom;a++) out.gradient[a]+=Gt[a];
                    for (int x=0;x<3;x++) for (int y=0;y<3;y++) V[x][y]+=Vt[x][y];
                }
            }
            if (firstEx) std::rethrow_exception(firstEx);
        }
        else
#endif
        for (const auto& [a,g] : work) point(a, g, out.gradient, V);
        out.virial=rmat3d_t(V[0][0],V[0][1],V[0][2], V[1][0],V[1][1],V[1][2], V[2][0],V[2][1],V[2][2]);
        return out;
    }
private:
    CellGeometry         itsCell;
    std::vector<rvec3_t> itsR;
    double               itsDPlane;
    int                  itsK;
    double               itsEps;
};

//! \brief An ATOM-CENTRED mesh must carry one site block per centre -- checked LIVE, not by \c assert.
//!
//! The Becke build's whole point is that each point belongs to a named centre and its weight already
//! carries that centre's partition \f$w_A\f$, which is what makes \f$\int w_A f\f$ a real atomic
//! integral (doc/OpenWork.md Step 0a).  A path that rebuilds the mesh and forgets to reopen the blocks
//! does not FAIL -- \c qcMesh::SiteIntegrals is simply never reached and the site-moment instrument
//! prints nothing.  That happened, in Release, on every imposed run, for as long as the invariant-mesh
//! filter walked orbits instead of sites; the \c assert that should have caught it is compiled out under
//! NDEBUG.  So this throws: a silently missing instrument is worse than a stopped run.
void RequireSiteBlocks(const qcMesh::Mesh& m, const char* who)
{
    if (m.NSites()>0) return;
    throw std::runtime_error(std::string(who)+": the atom-centred (Becke) mesh came out with NO site "
        "blocks.  Its per-centre partition is what makes an integrated site moment/charge possible, so a "
        "mesh without them silently disables that instrument -- some rebuild of the point list dropped "
        "MeshBuilder::BeginSite.");
}

//! \a angPerAtom (optional) = per-atom angular sets -- the §6a W2b SITE-ADAPTED path (rep atoms carry
//! site-group-invariant sets, partners the op-rotated copies; the ops-taking MakeMesh(UnitCell)
//! overload builds them).
//! Null = the historical shared single-orientation angular grid.
//! \a repOf (optional, with \a angPerAtom): per atom, its ORBIT REPRESENTATIVE (an atom index <= its own).  A member's angular
//! set is its representative's rotated by the op that maps one onto the other, point for point, and the Becke partition is
//! invariant under that space-group op -- so member point q carries EXACTLY the representative's weight and keep decision
//! (V1.22, anchor batch 2026-10-21).  Partitioned once per orbit: the mesh is orbit-consistent BY CONSTRUCTION (the eps-tail
//! drops used to be bit-different per member, so the caller post-filtered orbit-incomplete points) and ~|orbit| x cheaper.
qcMesh::Mesh MakePeriodicBeckeMesh(const UnitCell& cell, const qcMesh::MeshParams& mp,
                                   const std::vector<qcMesh::AngularMesh>* angPerAtom = nullptr,
                                   const std::vector<size_t>* repOf = nullptr)
{
    qchem::report::Timed timed("setup: becke mesh build");
    using qcMesh::BeckeCutoff;
    const int    k  =mp.beckeOrder;
    // Magnitude screen: the weights are eps-converged.  The tail bounds are worst-case-COLLINEAR (attained
    // only on the exact line through two images), so the TRUE neglected mass is far below eps -- while the
    // quadrature error itself sits at ~1e-4.  eps therefore has two decades of headroom, and it is the one
    // parameter that scales the dominant loop rather than shaving it: eps fixes |im|, the competitor count,
    // and the partition costs O(|P-set| x |im|) per point.
    //
    // 1e-6, MEASURED 2026-08-20 (doc/OpenWork.md Step 3).  This had only ever been probed TIGHTER -- 1e-10
    // was known to inflate the sets several-fold for no accuracy gain -- and the sweep in the other
    // direction is where the cost was:
    //     eps     images/pt (FCC)    BeckeCutoff calls    becke build (MnO, threaded)
    //     1e-8         1939               4.14e9                 36.95 s
    //     1e-7         1056               2.04e9                  ~18 s   (2.03x)
    //     1e-6          623               9.95e8                  8.33 s  (4.44x)
    //     1e-5          351               3.85e8                  --      (gate FAILS)
    // On MnO AFM-II VA the flip is invisible: Etot -52.48387019 -> -52.48387023 (9 s.f.), site moments
    // identical to 6 digits.  The BINDING gate is BeckeEquivalentSitesOwnEqualShares (site shares equal to
    // 1e-8 RELATIVE, a partition property and far sharper than any integration test); it holds at 1e-6 and
    // breaks at 1e-5, so the default keeps a full decade of margin against the gate that actually bites.
    // Unlike the Phi-table 4.6x this is a TOLERANCE trade, not a bit-identical restructuring -- the weights
    // do move, at ~1e-6 relative.  The value is MeshParams::beckeEps (default 1e-6); the environment override QCHEM_BECKE_EPS
    // (old name GPW_BECKE_EPS, deprecated) is applied where the facade builds the parameters (BeckeXCParams), not here.
    const double eps=mp.beckeEps;
    assert(k>=0);

    std::vector<rvec3_t> R;
    for (const Site& a : cell.Sites()) R.push_back(a.R);
    const size_t natom=R.size();
    assert(natom>0);

    // Minimum interplanar spacing d_i = 1/sqrt((M^-1)_ii): a point in cell n0 is at least s*dPlane from
    // anything at Chebyshev cell distance > s (s whole cell layers in between) -- the shell-tail floor.
    const Matrix3D<double>& A=cell.GetCellMatrix();
    Matrix3D<double> Minv=Invert(Transpose(A)*A);
    const double dPlane=1.0/sqrt(max(Minv(1,1),max(Minv(2,2),Minv(3,3))));

    qcMesh::RadialMesh  rad=qcMesh::MakeRadial(mp);   // one single-centre template reused (ShiftOrigin per atom)
    qcMesh::AngularMesh ang=qcMesh::MakeAngular(mp);

    static const bool census=qchem::Diagnostics::Enabled("becke_count");
    qcMesh::MeshBuilder out;
    // GPW_BECKE_COUNT totals over the whole mesh; each atom's slot-indexed census folds in below.
    struct { size_t pts=0, live=0, im=0, pset=0, retest=0, final_=0, bound=0,
                    dead=0, deadWork=0, imMax=0, psetMax=0, shellMax=0; } tot;
    assert(!repOf || (angPerAtom && repOf->size()==natom));
    std::vector<std::vector<char>> keepOf(natom);     // per REPRESENTATIVE: its keep decisions and kept weights (V1.22)
    std::vector<rvec_t>            kwtOf (natom);
    for (size_t ia=0; ia<natom; ia++)
    {
        // SITE BLOCK: this atom's surviving points, whose weights already carry ITS partition function
        // w_A(r) (kwt=wts*w below).  Recording the block is what lets a consumer ask the mesh for a real
        // ATOMIC integral -- Integral w_A(r) f(r) d3r = the block sum -- instead of sampling a field at a
        // guessed point and calling it a moment (doc/OpenWork.md Step 0a).  Costs one index per atom.
        out.BeginSite();
        qcMesh::Mesh am=qcMesh::ProductMesh(rad, angPerAtom ? (*angPerAtom)[ia] : ang);
        const rvec3vec_t u=am.Points();   // the TEMPLATE offsets, recorded per kept point (Mesh::OwnerOffsets: a pullback re-forms R_a + u exactly)
        am.ShiftOrigin(R[ia]);
        const rvec3vec_t& pts=am.Points();
        const rvec_t&     wts=am.Weights();
        // The per-point Becke partitions are INDEPENDENT (~1 ms each -- the eps-converged competitor
        // image series is the whole build cost, measured linear at any angular count), and each result
        // lands in its own indexed slot consumed in index order below -- so the threaded build is
        // BIT-IDENTICAL to the serial one at any thread count (no cross-point reductions, unlike the
        // GPW pair loops).  Hence threading it costs no anchor; QCHEM_OPENMP_THREADS (one knob
        // for every region, serial by default) sets the count -- purely resource control, never a numerics knob.
        const size_t nq=am.size();
        std::vector<char> keep(nq, 0);
        rvec3vec_t        kpt(nq);
        rvec_t            kwt(nq);
        std::vector<BeckePointCost> cost(census ? nq : 0);
        // A MEMBER of an orbit: the representative's decisions and weights, point for point; only the position is its own
        // (wrapped into the home cell exactly as the partition below wraps a kept point).
        if (repOf && (*repOf)[ia]!=ia)
        {
            const size_t r0=(*repOf)[ia];
            if (r0>ia || keepOf[r0].size()!=nq)
                throw std::logic_error("MakePeriodicBeckeMesh: an orbit member before its representative, or a member whose point "
                                       "count differs from its representative's");
            for (size_t q=0; q<nq; q++)
            {
                if (!keepOf[r0][q]) continue;
                const rvec3_t f=cell.ToFractional(pts[q]);
                out.Append(pts[q]-cell.ToCartesian(rvec3_t(floor(f.x),floor(f.y),floor(f.z))), kwtOf[r0][q], u[q]);
            }
            continue;
        }
        auto partition=[&](size_t q)
        {
            BeckePointCost pc;                       // recorded on every exit path when the census is on
            std::vector<qcMesh::BeckeImage> im;      // competitor images for this point (thread-local)
            long target=-1;
            const rvec3_t& r=pts[q];
            // THE ONE PARTITION (qcMesh::PeriodicBeckeWeight, lane B 2026-11-19): the build and its geometry derivative are one computation.
            const double w=qcMesh::PeriodicBeckeWeight(cell, R, dPlane, ia, r, k, eps, im, target, pc, census);
            if (census) { pc.im=unsigned(im.size()); cost[q]=pc; }
            if (w>0)   // wrapped into the home cell; slot-indexed so the threaded build stays ordered
            {
                const rvec3_t f=cell.ToFractional(r);
                keep[q]=1;
                kpt[q]=r-cell.ToCartesian(rvec3_t(int(floor(f.x)),int(floor(f.y)),int(floor(f.z))));
                kwt[q]=wts[q]*w;
            }
        };
#ifdef QCHEM_OPENMP
        const int nthreads=qchem::InnerThreads();     // THE one thread count (qchem.Parallel; 1 inside an outer-loop item)
        if (nthreads>1)
        {
            #pragma omp parallel for schedule(dynamic, 8) num_threads(nthreads)
            for (size_t q=0; q<nq; q++) partition(q);
        }
        else
            for (size_t q=0; q<nq; q++) partition(q);
#else
        for (size_t q=0; q<nq; q++) partition(q);
#endif
        for (size_t q=0; q<nq; q++) if (keep[q]) out.Append(kpt[q], kwt[q], u[q]);
        if (repOf) { keepOf[ia]=keep; kwtOf[ia]=kwt; }
        for (size_t q=0; q<cost.size(); q++)
        {
            const BeckePointCost& c=cost[q];
            tot.pts++;
            if (c.pset>0) tot.live++;             // reached the P-set double loop (not an early tail drop)
            tot.im+=c.im; tot.pset+=c.pset;
            tot.retest+=c.retest; tot.final_+=c.final_; tot.bound+=c.bound;
            tot.dead+=c.dead; tot.deadWork+=c.deadWork;
            tot.imMax=max(tot.imMax,size_t(c.im)); tot.psetMax=max(tot.psetMax,size_t(c.pset));
            tot.shellMax=max(tot.shellMax,size_t(c.shells));
        }
        // GPW_BECKE_ATOMS: the PER-ATOM breakdown.  Sum(w) is the partition's share of unity assigned to
        // this atom -- for two crystallographically EQUIVALENT atoms it must be identical, and so must the
        // kept/dropped counts.  The aggregate "[Becke grid] ... dropped N tail pts" line cannot show that,
        // which is why a site-dependent mesh defect could hide behind it (MnO 2026-08-11: the moment dies
        // on the atom at the cell CORNER and survives on the one at the cell CENTRE, and a rigid
        // translation of the whole crystal -- an exact symmetry -- moves the answer by 56 Ha).
        if (qchem::Diagnostics::Enabled("becke_atoms"))
        {
            report::FormatGuard keep_(std::cout);
            size_t nk=0; double wsum=0.0;
            for (size_t q=0; q<nq; q++) if (keep[q]) { ++nk; wsum+=kwt[q]; }
            std::cout<<"[Becke atom] ia="<<ia<<" R=("<<R[ia].x<<","<<R[ia].y<<","<<R[ia].z<<")"
                     <<" pts="<<nq<<" kept="<<nk<<" dropped="<<(nq-nk)
                     <<" Sum(w)="<<std::setprecision(10)<<wsum<<std::endl;
        }
    }
    qcMesh::Mesh mesh=out.take();
    mesh.InstallWeightResponse(std::make_shared<PeriodicBeckeResponse>(cell, R, dPlane, k, eps));   // its weights' derivative (the XC grid response)

    // Announce the grid -- the console line (the [GPW grid] family: one line, run start) is HOW a user
    // tells the Becke XC quadrature is in play at all; the report addendum lands in the `grids` section
    // next to the uniform ladder whenever a report run is open (inert otherwise, e.g. bare mesh tests).
    auto AngularName=[](qcMesh::AngularKind a)
    { switch (a) {
        case qcMesh::AngularKind::Lebedev:         return "Lebedev";
        case qcMesh::AngularKind::GaussLegendre: return "GaussLegendre";
        case qcMesh::AngularKind::None:          return "none (radial-only)"; }
      return "?"; };
    auto RadialName=[](qcMesh::RadialKind r)
    { switch (r) {
        case qcMesh::RadialKind::MHL:    return "MHL";
        case qcMesh::RadialKind::Log:    return "Log";
        case qcMesh::RadialKind::Linear: return "Linear";
        case qcMesh::RadialKind::DE2:    return "DE2"; }
      return "?"; };
    size_t nDirs=0;
    if (angPerAtom) for (const auto& a : *angPerAtom) nDirs=max(nDirs, a.size());
    else            nDirs=ang.size();
    size_t total=0;
    if (angPerAtom) for (const auto& a : *angPerAtom) total+=size_t(rad.size())*a.size();
    else            total=natom*size_t(rad.size())*ang.size();
    const size_t kept=mesh.size();
    const double angRot=angPerAtom ? 0.0 : mp.angRot;   // a site-adapted rule is the site's own: the rotation applies to the plain factory grid only
    std::cout<<"[Becke grid] natom="<<natom<<" radial="<<RadialName(mp.radial)<<" nR="<<mp.nRadial
             <<" alpha="<<(mp.radial==qcMesh::RadialKind::DE2 ? mp.de2_alpha : mp.mhl_alpha)<<" angular="<<(angPerAtom?"SITE-ADAPTED":AngularName(mp.angular))
             <<" L="<<mp.angularDegree<<" ("<<nDirs<<" dirs"<<(angPerAtom?" max/atom":"")<<")";
    if (angRot!=0.0) std::cout<<" angRot="<<angRot;
    std::cout<<" beckeOrder="<<k
             <<" points="<<kept<<" (dropped "<<total-kept<<" tail pts)"<<std::endl;
    if (census)
    {
        const double perPt=tot.live>0 ? double(tot.live) : 1.0;   // per LIVE point: the tail drops do no work
        const size_t calls=tot.retest+tot.final_+tot.bound;
        report::FormatGuard keep_(std::cout);
        std::cout<<std::setprecision(4)
                 <<"[Becke count] pts="<<tot.pts<<" live="<<tot.live
                 <<" images/pt="<<tot.im/perPt<<" (max "<<tot.imMax<<", shells<="<<tot.shellMax<<")"
                 <<" Pset/pt="<<tot.pset/perPt<<" (max "<<tot.psetMax<<")"
                 <<"\n[Becke count] BeckeCutoff/live-pt: final="<<tot.final_/perPt
                 <<" retest="<<tot.retest/perPt<<" bound="<<tot.bound/perPt
                 <<"  TOTAL="<<double(calls)
                 <<"  shares: final="<<100.0*tot.final_/double(max(calls,size_t(1)))<<"%"
                 <<" retest="<<100.0*tot.retest/double(max(calls,size_t(1)))<<"%"

                 <<"\n[Becke count] P-set members ending at P=0: "
                 <<100.0*tot.dead/double(max(tot.pset,size_t(1)))<<"% of members, "
                 <<100.0*tot.deadWork/double(max(tot.final_,size_t(1)))<<"% OF THE FINAL-LOOP WORK => "
                 <<1.0/(1.0-double(tot.deadWork)/double(max(tot.final_,size_t(1))))<<"x if avoided"<<std::endl;
    }
    // ONE LEVEL IN (user ruling 2026-11-02): the CALLER opened a report::Scope naming itself; the builder writes its fields
    // there and picks no key.  Was a root-anchored EmitAt("grids","becke") -- one key, so a second Becke mesh overwrote it.
    namespace rpt=qchem::report;
    rpt::Set("kind", "Becke"); rpt::Set("natom", (long)natom); rpt::Set("radial", RadialName(mp.radial));
    rpt::Set("nRadial", mp.nRadial); if (mp.radial==qcMesh::RadialKind::DE2) { rpt::Set("de2_alpha", mp.de2_alpha); rpt::Set("de2_rMax", mp.de2_rMax); } else rpt::Set("mhl_alpha", mp.mhl_alpha);
    rpt::Set("angular", angPerAtom ? "SiteAdapted" : AngularName(mp.angular)); rpt::Set("L", mp.angularDegree); rpt::Set("angRot", angRot);
    rpt::Set("nDirs", (long)nDirs); rpt::Set("beckeOrder", k); rpt::Set("points", (long)kept);
    rpt::Set("dropped", (long)(total-kept)); rpt::Set("eps", eps);
    return mesh;
}

} //anon

// §6a W2b: the SITE-ADAPTED periodic Becke mesh.  Per atom ORBIT of the imposed ops, the representative
// atom carries a site-group-INVARIANT angular set (MakeInvariantAngularMesh under its CARTESIAN stabilizer
// A·W·A⁻¹) and every symmetry partner carries the op-rotated copy (the atom-fold's edge op) -- so the whole
// mesh maps onto itself under every op BY CONSTRUCTION (the T2 precondition), with the fuzzy-Voronoi
// partition computed on the final point set (a geometric function of the atom distances, so the weights
// inherit the invariance).
//! \brief GPW_MESH_ORTHO=M: HOW ORTHOGONAL ARE PLANE WAVES UNDER THIS MESH'S QUADRATURE?
//!
//! WHY IT IS ONE FUNCTION AND NOT A MATRIX.  The discrete overlap of two plane waves on a mesh,
//! \f$S_{GG'}=\sum_g w_g e^{i(G'-G)\cdot r_g}\f$, depends only on the DIFFERENCE, so the whole
//! non-orthogonality is one scalar function \f$T(\Delta G)\f$ over the difference set -- O(npts) per
//! \f$\Delta G\f$, with no \f$N_G^2\f$ object anywhere.  \f$T(\Delta G)/\Omega\f$ IS the mesh's
//! quadrature error for a plane wave of wavevector \f$\Delta G\f$, and \f$T(0)/\Omega=1\f$ exactly by
//! construction -- a free correctness check on the census itself.
//!
//! WHAT IT IS FOR.  A plane-wave \c VxcFit on a Becke grid is asserted out (item I3 in Hamiltonians.C) and
//! the projection \f$\langle c_G|v_{xc}\rangle=\sum_g w_g v_{xc}(r_g)e^{-iG\cdot r_g}\f$ is a perfectly
//! ordinary real-space quadrature that takes any mesh (user, 2026-08-22) -- but on a UNIFORM raster the
//! plane waves are DISCRETELY ORTHOGONAL, so the fit is a projection with no metric, while on any other
//! mesh it becomes a genuine fit with an overlap to invert.  This measures how far from orthogonal.
//!
//! THE EXPECTATION, so the measurement can refute it: a Becke mesh is built for the OPPOSITE integrand --
//! dense radial shells at the nuclei, sparse in the interstitial, tuned for sharp atom-centred cusps.  A
//! plane wave has uniform amplitude everywhere and oscillates fastest exactly where the mesh is thinnest,
//! and through a radial shell it is an ANGULAR oscillation -- so the error should be ~0 at small
//! \f$|\Delta G|\f$ (a near-constant integrand) and degrade once \f$|\Delta G|\f$ outruns the Lebedev
//! order.  Reported binned in \f$|\Delta G|\f$ so the onset is visible rather than averaged away.
void ReportMeshPlaneWaveOrtho(const UnitCell& cell, const qcMesh::Mesh& mesh)
{
    const auto e=qchem::Diagnostics::Value("mesh_ortho");
    if (!e) return;
    const int M=std::max(1, std::atoi(e->c_str()));            // half-width of the integer Delta-G box (default 1 -> set 4)
    const rvec3vec_t& R=mesh.Points();
    const rvec_t&     W=mesh.Weights();
    if (R.size()==0) return;
    double Omega=0.0; for (size_t g=0; g<W.size(); ++g) Omega+=W[g];
    const CellGeometry B=cell.MakeReciprocal();

    // Bin in |Delta G| (a.u.^-1).  Geometric-ish edges: the interesting behaviour is an ONSET, not a slope.
    constexpr size_t kNBin=8;
    const double edge[kNBin]={0.5,1.0,1.5,2.0,3.0,4.0,6.0,1e300};
    double binMax[kNBin]={0,0,0,0,0,0,0,0}, binSum[kNBin]={0,0,0,0,0,0,0,0};
    size_t binN[kNBin]={0,0,0,0,0,0,0,0};
    double t0=0.0;
    for (int mx=-M; mx<=M; ++mx) for (int my=-M; my<=M; ++my) for (int mz=-M; mz<=M; ++mz)
    {
        const rvec3_t G=B.ToCartesian(rvec3_t(double(mx),double(my),double(mz)));
        const double  g2=G*G;
        double re=0.0, im=0.0;
        for (size_t q=0; q<R.size(); ++q)
        {
            const double ph=G*R[q];
            re+=W[q]*cos(ph); im+=W[q]*sin(ph);
        }
        const double t=sqrt(re*re+im*im)/Omega;
        if (g2<1e-20) { t0=t; continue; }             // Delta G = 0: the T(0)/Omega == 1 self-check
        const double gg=sqrt(g2);
        size_t b=0; while (b+1<kNBin && gg>=edge[b]) b++;
        binMax[b]=max(binMax[b],t); binSum[b]+=t; binN[b]++;
    }
    std::cout<<"[mesh ortho] npts="<<R.size()<<"  T(0)/Omega="<<t0
             <<" (must be 1)  |Delta G| box M="<<M<<std::endl;
    const char* lbl[kNBin]={"<0.5","0.5-1","1-1.5","1.5-2","2-3","3-4","4-6",">6"};
    std::cout<<"[mesh ortho] |T(dG)|/Omega  (0 = orthogonal):";
    for (size_t b=0;b<kNBin;b++)
        if (binN[b]) std::cout<<"  "<<lbl[b]<<": max="<<binMax[b]<<" mean="<<binSum[b]/double(binN[b])
                              <<" (n="<<binN[b]<<")";
    std::cout<<std::endl;
}

// WHAT A CELL'S MESH FORCES, from the one descriptor (qcMesh::Forces; doc/Pins.md pin 30, Group face step 3b) in the cell's frame: its
// lattice, and no centre -- a cell has a lattice of sites, and no field periodic on a lattice is spherical about one point (only a
// constant is both).  So the radial-only set (angular none) cannot exist here: a CONSTRUCTION REFUSAL, the descriptor's own reason.
// \a mp is the mesh AS BUILT: an Auto kind reaching the uniform branch is stamped as the uniform raster it became.
static Symmetry::GroupPtr CellMeshForces(const UnitCell& cell, const qcMesh::MeshParams& mp)
{
    const auto H=qcMesh::Forces(mp, Symmetry::Frame{.lattice=cell.GetCellMatrix(), .centre={}});
    if (!H) throw std::invalid_argument("the periodic cell mesh: "+H.Error()+" (doc/Pins.md pin 30)");
    return H.Value();
}
static qcMesh::MeshParams AsUniform(qcMesh::MeshParams mp) {mp.cellKind=qcMesh::CellKind::RecipeStated(qcMesh::UnitCellKind::Uniform); return mp;}

qcMesh::Mesh MakeMesh(const UnitCell& cell, const qcMesh::MeshParams& mp, const std::vector<Symmetry::SymOp>& ops)
{
    (void)CellMeshForces(cell, mp);   // refuse what cannot exist on a cell, before any point is built
    namespace SL = Symmetry::Lattice_3D;
    assert(!ops.empty());

    // UNIFORM kind: there is no site-adapted construction to do -- build the generic mesh and group-average
    // it onto the op orbits (W1; a no-op dedup when the grid is already commensurate).  This branch used to
    // sit in the CALLER (GPW_IBS's XC-quadrature factory), which meant a BASIS was switching on cellKind to
    // decide how a STRUCTURE builds its own mesh.  The structure owns the mesh; it owns the strategy.
    if (mp.cellKind!=qcMesh::UnitCellKind::Becke)
    {
        qcMesh::Mesh m=MakeInvariant(MakeMesh(cell, mp), cell.GetCellMatrix(), ops, kMeshMatchTol);
        m.InstallForces(CellMeshForces(cell, AsUniform(mp)));   // the orbit-averaged raster is still this cell's lattice of points
        ReportMeshPlaneWaveOrtho(cell, m);
        return m;
    }

    // BECKE kind: site-adapted, invariant by construction (W2b) -- the rest of this function.
    // D5 (anchor batch, 2026-10-21): the ANGULAR request, resolved and said.  The site-adapted builder consumes only the
    // DEGREE: a stated scheme (Lebedev / GaussLegendre) cannot be honoured under an imposed group and is REPLACED by the
    // invariant rule of the same degree -- say so once, as a warning (it used to be silent).  An UNSET degree (<= 0) takes the
    // calibrated Becke default.
    const int L = mp.angularDegree>0 ? mp.angularDegree : 29;   // 29 = qcMesh::BeckeXCParams' calibrated default (qcMesh XCPolicy, not imported here)
    {
        static bool warned=false;
        if (!warned)
        {
            warned=true;
            std::cout<<"[Becke grid] ⚠ imposed symmetry: the requested "
                     <<(mp.angular==qcMesh::AngularKind::Lebedev ? "Lebedev" : "GaussLegendre")<<" angular rule"
                     <<(mp.angularDegree>0 ? "" : " (degree unset -> the Becke default)")
                     <<" is REPLACED by the site-adapted invariant rule of degree "<<L<<std::endl;
        }
    }

    std::vector<rvec3_t> F;                         // fractional atom sites
    for (const Site& a : cell.Sites()) F.push_back(cell.ToFractional(a.R));
    SL::Fold af = SL::FoldPointsPeriodic(F, ops, 1e-6);   // atom orbits + the rep->partner edge ops

    const Matrix3D<double>& A = cell.GetCellMatrix();
    const Matrix3D<double> Ainv = Invert(A);
    auto cart  = [&](const Matrix3D<double>& W) {return A*W*Ainv;};   // fractional op -> Cartesian rotation
    auto fixes = [&](const SL::SymOp& op, const rvec3_t& f)           // torus site-stabilizer test
    {
        rvec3_t im = op.W*f + op.tau;
        double dx = im.x-f.x, dy = im.y-f.y, dz = im.z-f.z;
        dx -= floor(dx+0.5); dy -= floor(dy+0.5); dz -= floor(dz+0.5);
        return dx*dx+dy*dy+dz*dz < 1e-12;
    };

    // The mixed rule's ATOM-SPECIFIC bad-direction list (§6a W2): the atom's actual bond directions
    // -- unit vectors to its nearest-neighbour images (first shell, 5% width; the +-1 cell range
    // covers the nearest neighbour of any reasonable primitive cell).  Only these are poison for a
    // quadrature direction; special orbits pointing BETWEEN bonds stay admissible.
    auto bondDirs = [&](const rvec3_t& fr)
    {
        const rvec3_t r0 = cell.ToCartesian(fr);
        std::vector<std::pair<double,rvec3_t>> nb;
        double dmin = 1e300;
        for (const Site& a : cell.Sites())
            for (int i = -1; i <= 1; ++i)
                for (int j = -1; j <= 1; ++j)
                    for (int k = -1; k <= 1; ++k)
                    {
                        rvec3_t dR = a.R + cell.ToCartesian(rvec3_t(i,j,k)) - r0;
                        double  d  = norm(dR);
                        if (d < 1e-6) continue;          // the atom itself
                        nb.push_back({d, dR/d});
                        if (d < dmin) dmin = d;
                    }
        std::vector<rvec3_t> dirs;
        for (const auto& [d, u] : nb) if (d < dmin*1.05) dirs.push_back(u);
        return dirs;
    };

    std::vector<qcMesh::AngularMesh> ang(F.size());
    // Priced with the rest of the imposed build (1.1(b)); stopped explicitly just before the Becke build,
    // which has its own bucket.
    auto angTimer = std::make_unique<qchem::report::Timed>("setup: site-adapted angular sets (W2b)");
    for (size_t o = 0; o < af.Reps(); ++o)
    {
        const rvec3_t& fr = F[af.repRaw[o]];
        std::vector<Matrix3D<double>> siteOps;
        for (const auto& op : ops) if (fixes(op, fr)) siteOps.push_back(cart(op.W));
        qcMesh::Mesh aq = MakeInvariantAngularMesh(siteOps, L, bondDirs(fr));
        for (auto [mi, oi] : af.members[o])
        {
            assert(oi>=0 && "site-adapted Becke mesh: the imposed op set must be a group (edge op per partner)");
            const Matrix3D<double> Rrot = cart(ops[oi].W);
            rvec3vec_t d(aq.size());
            rvec_t     w(aq.size());
            for (size_t i = 0; i < aq.size(); ++i) {d[i] = Rrot*aq.Points()[i]; w[i] = aq.Weights()[i];}
            ang[mi] = qcMesh::AngularMesh(std::move(d), std::move(w));
        }
    }
    angTimer.reset();                       // the site-adapted angular sets are built; stop that bucket
    // V1.22: one partition per ATOM ORBIT, copied to its members point for point (see MakePeriodicBeckeMesh) -- the mesh is
    // orbit-consistent by construction.  This replaced the post-filter that dropped orbit-INCOMPLETE eps-borderline points
    // (the builder's per-point tail drops were bit-different across an orbit) with a SECOND full orbit fold over ~97k points.
    std::vector<size_t> repOf(F.size());
    for (size_t o = 0; o < af.Reps(); ++o)
        for (auto [mi, oi] : af.members[o]) repOf[mi] = size_t(af.repRaw[o]);
    qcMesh::Mesh mb = MakePeriodicBeckeMesh(cell, mp, &ang, &repOf);
    RequireSiteBlocks(mb, "MakeMesh(UnitCell, imposed, Becke)");
    ReportMeshPlaneWaveOrtho(cell, mb);   // the BECKE arm: built for cusps, asked about plane waves
    return mb;
}

// A periodic cell's real-space integration mesh.  Two kinds (mp.cellKind):
// Uniform -- a UNIFORM grid over the cell: n=mp.nUniform points per axis at cell-fractional midpoints
// f = ((i+1/2)/n, (j+1/2)/n, (k+1/2)/n), mapped to Cartesian by r = A f (ToCartesian), each carrying equal
// weight Omega/n^3.  For a cell-periodic integrand the midpoint rule is the natural (and, for smooth
// periodic fields, exponentially convergent) quadrature.  This is the working lattice mesh that lets the
// geometry-neutral PP terms (PP_Local/PP_NonLocal) run on a real-space lattice basis.  (Plane-wave DFT
// integrates in G-space on the basis's own grid, so it never asks for this.)
// Becke -- the atom-centred periodic fuzzy-Voronoi quadrature above (dense radial near each nucleus,
// cheap diffuse tails): the XC quadrature for pointwise-nonlinear, sharp-at-the-core fields.
qcMesh::Mesh MakeMesh(const UnitCell& cell, const qcMesh::MeshParams& mp)
{
    (void)CellMeshForces(cell, mp);   // refuse what cannot exist on a cell, before any point is built
    if (mp.cellKind==qcMesh::UnitCellKind::Becke)
    {
        qcMesh::Mesh mB=MakePeriodicBeckeMesh(cell,mp);
        RequireSiteBlocks(mB, "MakeMesh(UnitCell, free, Becke)");
        ReportMeshPlaneWaveOrtho(cell, mB);   // the FREE (non-imposed) Becke arm
        return mB;
    }

    // Points per axis: DERIVE from a physical grid cutoff when mp.eCut>0 (the GPW / Nyquist path), else the
    // manual mp.nUniform.  The Nyquist relation itself is qcMesh::UniformDivisions -- named once beside the
    // two MeshParams fields it relates, because the resolution CHECKS (is this mesh fine enough for the
    // basis? which grid is cheaper?) need the same mapping and must not re-derive it.
    const int n = mp.eCut>0.0 ? qcMesh::UniformDivisions(cell.SizingGeometry().GetMaximumCellEdge(), mp.eCut) : mp.nUniform;   // the frozen rule
    assert(n>0);
    const double w=cell.GetCellVolume()/(double(n)*n*n);   // equal weight Omega/n^3
    rvec3vec_t R(size_t(n)*n*n);
    rvec_t     W(size_t(n)*n*n);
    size_t q=0;
    for (int i=0; i<n; i++)
        for (int j=0; j<n; j++)
            for (int k=0; k<n; k++,q++)
            {
                R[q]=cell.ToCartesian(rvec3_t((i+0.5)/n, (j+0.5)/n, (k+0.5)/n));   // r = A f, fractional midpoints
                W[q]=w;
            }

    // Announce, like the Becke branch -- so a run's report always says WHICH cell quadrature is in play -- and LIKE it,
    // ONE LEVEL IN (user ruling 2026-11-02): the caller (PP_Local, PP_NonLocal, the GPW PP quadrature, the XC quadrature)
    // opened a report::Scope naming itself, and these fields land there.  Was a root-anchored grids.unitCellUniform.
    cell.SetUniformGridFields(vec3_t<int>(n,n,n), mp.eCut);
    qcMesh::Mesh m(std::move(R), std::move(W));
    // The points ARE the cell's midpoint raster in (i,j,k) row-major order, so a Fourier series on this cell is sampled
    // here by the exact inverse FFT, not the NUFFT default (the rho-tilde sampler, OpenWork §4; user ruling 2026-11-04).
    m.InstallFourierSampling(std::make_shared<const MidpointRasterSampling>(ivec3_t(n,n,n),
        cell.ToCartesian(rvec3_t(1,0,0)), cell.ToCartesian(rvec3_t(0,1,0)), cell.ToCartesian(rvec3_t(0,0,1))));
    m.InstallForces(CellMeshForces(cell, AsUniform(mp)));   // a field on the raster is periodic on this lattice by construction (pin 30)
    ReportMeshPlaneWaveOrtho(cell, m);   // the UNIFORM control: exact orthogonality up to Nyquist
    return m;
}

// ---- The Becke molecular integration mesh (module-internal; MakeMesh(Atom) and MakeMesh(Molecule)) --------
// A. D. Becke, J. Chem. Phys. 88(4), 2547 (1988).  Cell function for atom i at point r:
//   p_i(r) = prod_{j != i} s(mu_ij),   mu_ij = (|r-R_i| - |r-R_j|) / |R_i - R_j|,   s(mu) = 1/2(1 - f_k(mu)).
// mu is bounded in [-1,1] for ANY nonzero separation (reverse triangle inequality); the only singularity is
// exactly coincident atoms (R_ij = 0), where we set mu = 0 (so s = 1/2 -> a coincident dimer integrates to
// the single-atom result).  Serves both Atom (natom==1: just the shifted product grid) and Molecule.
// The smoothing polynomial s(mu) is the shared qcMesh::BeckeCutoff (also used by the periodic
// UnitCell Becke partition).
using qcMesh::BeckeCutoff;

// A STATED uniform cell mesh is refused: it would impose a translation group this structure lacks (a representation may not
// impose a symmetry the structure lacks -- UTHamiltonian VxcProductSpace).  \c Auto and \c Becke both mean the atom-centred mesh.
// Until 2026-10-31 \c cellKind was ignored here, so a uniform request ran Becke silently.
static void RefuseCellMesh(const qcMesh::MeshParams& mp, const char* what)
{
    if (mp.cellKind==qcMesh::UnitCellKind::Uniform)
        throw std::invalid_argument(std::string(what)+": a uniform CELL mesh imposes a translation group; this structure has none "
                                    "(point-group symmetry only) -- To fix this, please state cellKind Becke or Auto -- the atom-centred mesh.");
}

// THE MESH REPORTS ITSELF as it is built (OpenWork §3 "A molecule's mesh never reports itself", 2026-11-02), ONE LEVEL IN:
// the CALLER opens a report::Scope naming itself (the XC quadrature, a PP term, the auxiliary fit, the residual measure ...)
// and this builder writes cursor-relative fields there -- the periodic builder's field names (UnitCell.C), never a key or a
// path of its own (user ruling: the caller sets the context, the callee reports one level in).  Inert outside a run.
static const char* RadialName(qcMesh::RadialKind r)
{
    switch (r)
    {
        case qcMesh::RadialKind::MHL:    return "MHL";
        case qcMesh::RadialKind::Log:    return "Log";
        case qcMesh::RadialKind::Linear: return "Linear";
        case qcMesh::RadialKind::DE2:    return "DE2";
    }
    return "?";
}
static const char* AngularName(qcMesh::AngularKind a)
{
    switch (a)
    {
        case qcMesh::AngularKind::Lebedev:       return "Lebedev";
        case qcMesh::AngularKind::GaussLegendre: return "GaussLegendre";
        case qcMesh::AngularKind::None:          return "none (radial-only)";
    }
    return "?";
}
static void ReportMolecularMesh(size_t natom, const qcMesh::MeshParams& mp, size_t nDirs, size_t total, size_t kept)
{
    namespace rpt=qchem::report;
    rpt::Set("kind", "Becke"); rpt::Set("natom", (long)natom); rpt::Set("radial", RadialName(mp.radial));
    rpt::Set("nRadial", mp.nRadial); if (mp.radial==qcMesh::RadialKind::DE2) { rpt::Set("de2_alpha", mp.de2_alpha); rpt::Set("de2_rMax", mp.de2_rMax); } else rpt::Set("mhl_alpha", mp.mhl_alpha); rpt::Set("angular", AngularName(mp.angular));
    rpt::Set("L", mp.angularDegree); rpt::Set("nDirs", (long)nDirs); rpt::Set("beckeOrder", mp.beckeOrder);
    rpt::Set("points", (long)kept); rpt::Set("dropped", (long)(total-kept));
}

//! \brief THE MOLECULAR BECKE MESH's WEIGHT RESPONSE (qcMesh::MeshWeightResponse; lane B 2026-11-22, the molecular twin of the periodic
//! one in UnitCell.C).  Installed by AtomCentredMesh (MakeMesh for an Atom or a Molecule), which alone knows how the weights were made: each kept point g of site a carries
//! w_g = w_q W_a(R_a + u_g), w_q its template weight, W_a the Becke partition (qcMesh::BeckePartitionWeight, order k).  For f held per point,
//! d(Σ_g f_g w_g)/dR_D = Σ_g f_g w_q dW_a/dR_D, the point riding its owner (qcMesh::BeckePartitionGradient, which sums to zero over D), and
//! the virial Σ_D (d/dR_D) ⊗ R_D.  ONE home for the partition: the build and this derivative call the same leaf.  A consumer is any term
//! whose energy is Σ_g w_g f_g on this mesh at fixed f: the δ XC energy (f = e_xc), and the outer quadrature of a δ Hartree.
class MolecularBeckeResponse : public qcMesh::MeshWeightResponse
{
public:
    MolecularBeckeResponse(std::vector<rvec3_t> R, int k) : itsR(std::move(R)), itsK(k) {}
    Result Pullback(const qcMesh::Mesh& m, const rvec_t& f) const override;
private:
    std::vector<rvec3_t> itsR;
    int itsK;
};

MolecularBeckeResponse::Result MolecularBeckeResponse::Pullback(const qcMesh::Mesh& m, const rvec_t& f) const
{
    const size_t natom=itsR.size();
    if (f.size()!=m.size()) throw std::invalid_argument("MolecularBeckeResponse: the field is "+std::to_string(f.size())+" for "+std::to_string(m.size())+" points");
    if (m.OwnerOffsets().size()!=m.size() || m.NSites()!=natom)
        throw std::logic_error("MolecularBeckeResponse: the mesh carries no owner offsets or not one site block per atom");
    qchem::report::Timed timed("force: becke weight response");
    Result out{rvec3vec_t(natom, rvec3_t(0,0,0)), rmat3d_t(0,0,0, 0,0,0, 0,0,0)};
    if (natom==1) return out;   // one atom: W = 1 everywhere, no response
    const rvec3vec_t& U=m.OwnerOffsets();
    const rvec_t& W=m.Weights();
    for (size_t a=0; a<natom; a++)
        for (size_t g=m.SiteBegin(a); g<m.SiteEnd(a); g++)
        {
            if (f[g]==0.0) continue;
            const rvec3_t r=itsR[a]+U[g];
            const double wP=qcMesh::BeckePartitionWeight(itsR, a, r, itsK);
            // THE PLATEAUS (as the periodic response; measured there 2026-11-19): where W is 1 or 0 to 1e-14 its derivative is of that order.
            if (wP>1.0-1e-14 || wP<1e-14) continue;
            const double s=f[g]*W[g]/wP;   // f_g w_q
            const std::vector<rvec3_t> dW=qcMesh::BeckePartitionGradient(itsR, a, r, itsK);
            for (size_t D=0; D<natom; D++) out.gradient[D]+=s*dW[D];
        }
    for (size_t D=0; D<natom; D++)
    {
        const double gg[3]={out.gradient[D].x, out.gradient[D].y, out.gradient[D].z}, xx[3]={itsR[D].x, itsR[D].y, itsR[D].z};
        for (int x=0;x<3;x++) for (int y=0;y<3;y++) out.virial(x+1,y+1)+=gg[x]*xx[y];
    }
    return out;
}

//! \brief THE POINT-GROUP BUILD of the atom-centred mesh (V-FOLD-NEUTRAL step 1): what AtomCentredMesh needs to build a mesh invariant
//! under a molecule's ops.  \c frame rotates the ONE shared angular set (its columns are the frame axes in Cartesian coordinates);
//! \c repOf names each atom's orbit representative (an atom index <= its own) and \c edge the Cartesian op that maps the representative
//! onto it.  A member's points are its representative's, carried by that op POINT FOR POINT with the weight copied: the Becke partition
//! depends only on distances, so the op leaves it invariant, and copying makes the orbit's weights agree bit for bit (the solid's V1.22
//! lesson, where an independently built member dropped eps-borderline points differently).
struct MolecularOrbitBuild
{
    Matrix3D<double>              frame;
    std::vector<size_t>           repOf;
    std::vector<Matrix3D<double>> edge;
};

static qcMesh::Mesh AtomCentredMesh(const Structure& atoms, const qcMesh::MeshParams& mp, const MolecularOrbitBuild* sym=nullptr)
{
    RefuseCellMesh(mp, "the atom-centred (molecule / atom) mesh");
    // WHAT THIS MESH FORCES, from the one descriptor (qcMesh::Forces; doc/Pins.md pin 30, Group face step 3b) in the structure's own frame
    // (its group's: an atom offers its site, a molecule's finite point group offers no centre).  Only the radial-only set (angular none)
    // forces anything, O(3) about the one centre, so it exists only on the atom: elsewhere a CONSTRUCTION REFUSAL, the descriptor's reason.
    Symmetry::GroupPtr forced=Symmetry::TrivialGroup();
    if (mp.angular==qcMesh::AngularKind::None)
    {
        const auto H=qcMesh::Forces(mp, atoms.MakeSymmetry(Symmetry::kAdmissionTolerance)->MakeFrame());
        if (!H) throw std::invalid_argument("the atom-centred mesh: "+H.Error()+" (doc/Pins.md pin 30)");
        forced=H.Value();
    }
    const int beckeOrder=mp.beckeOrder;
    assert(beckeOrder>=0);
    size_t natom=atoms.GetNumAtoms();

    std::vector<rvec3_t> R(natom);
    for (size_t i=0; i<natom; i++) R[i]=atoms.Sites()[i].R;

    qcMesh::RadialMesh  rad=qcMesh::MakeRadial(mp);   // one single-center template reused (ShiftOrigin per atom)
    qcMesh::AngularMesh ang=qcMesh::MakeAngular(mp);
    if (sym)   // the shared set, rotated into the point group's frame: invariant under every op of D2h type there (Lebedev is Oh-invariant)
    {
        rvec3vec_t d(ang.size());
        for (size_t i=0; i<ang.size(); i++) d[i]=sym->frame*ang.Dirs()[i];
        ang=qcMesh::AngularMesh(std::move(d), ang.W());
    }

    // The representatives' kept points (owner offset u, weight), for their orbit members to copy: filled only on the symmetric build.
    std::vector<std::vector<std::pair<rvec3_t,double>>> kept(sym ? natom : 0);
    qcMesh::MeshBuilder out;
    for (size_t ia=0; ia<natom; ia++)
    {
        out.BeginSite();   // one block per atom, in atom order (SiteIntegrals, the radial Poisson operator)
        if (sym && sym->repOf[ia]!=ia)   // an orbit member: its representative's points, carried by the edge op (the rep comes first)
        {
            const Matrix3D<double>& W=sym->edge[ia];
            for (const auto& [u, w] : kept[sym->repOf[ia]])
            {
                const rvec3_t uu=W*u;
                out.Append(R[ia]+uu, w, uu);
            }
            continue;
        }
        qcMesh::Mesh am=qcMesh::ProductMesh(rad,ang);
        am.ShiftOrigin(R[ia]);
        const rvec3vec_t& pts=am.Points();
        const rvec_t&     wts=am.Weights();
        for (size_t q=0; q<am.size(); q++)
        {
            const rvec3_t& r=pts[q];
            // The symmetric build drops ZERO-WEIGHT points: a radial node at r = 0 stacks one per direction on the nucleus, and coincident
            // points make the orbit fold ambiguous (one point would join two orbits).  They integrate nothing, so nothing else changes.
            if (sym && wts[q]==0.0) continue;
            if (natom==1) { out.Append(r, wts[q], r-R[ia]); continue; }

            // ONE home for the partition weight: the force's grid response differentiates exactly this (qchem.Mesh.BeckePartition)
            const double w=qcMesh::BeckePartitionWeight(R, ia, r, beckeOrder);
            if (w<=0.0) continue;
            out.Append(r, wts[q]*w, r-R[ia]);   // the owner offset: the point rides its atom (the weight response)
            if (sym) kept[ia].push_back({r-R[ia], wts[q]*w});
        }
    }
    qcMesh::Mesh mesh=out.take();
    mesh.InstallWeightResponse(std::make_shared<MolecularBeckeResponse>(R, beckeOrder));   // its weights' derivative (the XC grid response)
    mesh.InstallForces(forced);
    ReportMolecularMesh(natom, mp, ang.size(), rad.size()*ang.size()*natom, mesh.size());
    return mesh;
}

// ---- The atom-centred overloads -----------------------------------------------------------------------------------------------------
// Every atom-centred mesh passes through the factory, so the per-site block contract is checked on EVERY one of them (the plan's §7
// gain): until step 2 RequireSiteBlocks ran on the two periodic Becke paths only, and a molecule or atom mesh that lost
// MeshBuilder::BeginSite would have silently disabled the site-moment instrument.  O(sites); a throw that survives NDEBUG.
qcMesh::Mesh MakeMesh(const Atom& atom, const qcMesh::MeshParams& mp)
{
    qcMesh::Mesh m=AtomCentredMesh(atom, mp);
    RequireSiteBlocks(m, "MakeMesh(Atom)");
    return m;
}

qcMesh::Mesh MakeMesh(const Molecule& mol, const qcMesh::MeshParams& mp)
{
    qcMesh::Mesh m=AtomCentredMesh(mol, mp);
    RequireSiteBlocks(m, "MakeMesh(Molecule)");
    return m;
}

// ------------------------------------------------------------------ the symmetric overload's hooks (StructureMeshMaker<S>)

// The orbit fold, moved from CrystalImposition::MakeXCMesh (controller + user rulings 2026-11-22): the geometry it needs -- the cell
// matrix, the flip-fixed test's fractional coordinates -- is the cell's, so it is the maker's, never the source's.
qcMesh::FoldedMesh FoldOrbits(const UnitCell& cell, qcMesh::Mesh mesh, const std::vector<Symmetry::SymOp>& ops)
{
    const double tol = qchem::kMeshMatchTol;   // the ONE mesh-point coincidence tolerance (R2.12)
    const Matrix3D<double>& A = cell.GetCellMatrix();
    // ★ 1.1(b): the XC quadrature build is 23.75 s/call EXCLUSIVE of the mesh build and the Φ tables, and it is the largest
    // non-threading block in a GPW run.  These two buckets are the only candidates left inside it -- both are point-matching
    // passes over ~97k mesh points against 24 ops.
    auto fold = [&]{ qchem::report::Timed timed("setup: XC mesh orbit fold (FoldMesh)");
                     return FoldMesh(mesh, A, ops, tol); }();
    std::cout << "[Becke grid] imposed symmetry (" << ops.size() << " ops): invariant mesh "
              << mesh.size() << " points in " << fold.Reps()
              << " orbits; rho star-averaged each iteration (doc/SymmetryUpgradePlan.md 6a)" << std::endl;
    const bool magnetic=std::any_of(ops.begin(), ops.end(), [](const Symmetry::SymOp& op){return op.sigma!=Symmetry::SpinAction::None;});
    if (!magnetic)
        return {std::make_shared<const qcMesh::Mesh>(std::move(mesh)), std::move(fold)};

    // The σ tags (parallel to ops -- the fold's edge opIndex indexes them) and the odd-field zero flags: FRACTIONAL mesh
    // coordinates for the torus-metric fixed-point test.
    std::vector<Symmetry::SpinAction> sigmas;
    sigmas.reserve(ops.size());
    for (const auto& op : ops) sigmas.push_back(op.sigma);
    auto flags = [&]{
        qchem::report::Timed timed("setup: XC mesh Shubnikov flip-fixed test");
        std::vector<rvec3_t> frac;
        frac.reserve(mesh.size());
        const Matrix3D<double> Ainv = Invert(A);
        for (size_t i=0;i<mesh.size();++i) frac.push_back(Ainv*mesh.Points()[i]);
        return Symmetry::Lattice_3D::FlipFixedPointsPeriodic(frac, ops, tol); }();
    size_t nFixed=0; for (char c : flags) if (c) nFixed++;
    std::cout << "[Becke grid] Shubnikov: " << nFixed << " flip-fixed mesh points (m==0 there exactly)"
              << std::endl;
    return {std::make_shared<const qcMesh::Mesh>(std::move(mesh)), std::move(fold),
            std::move(sigmas), std::move(flags)};
}

qcMesh::Mesh MakeMesh(const Atom& atom, const qcMesh::MeshParams& mp, const std::vector<Symmetry::SymOp>& ops)
{
    if (mp.angular!=qcMesh::AngularKind::None)
        throw std::invalid_argument("MakeMesh(Atom, "+std::to_string(ops.size())+" ops): an angular mesh on an atom is not invariant under "
                                    "the atom's rotations (a finite direction set claims a symmetry it does not carry) -- To fix this, please "
                                    "state the radial-only set (xcMesh angular None), whose angular part is analytic");
    return MakeMesh(atom, mp);
}
qcMesh::FoldedMesh FoldOrbits(const Atom&, qcMesh::Mesh mesh, const std::vector<Symmetry::SymOp>&)
{
    return {std::make_shared<const qcMesh::Mesh>(std::move(mesh)), {}};   // one direction per shell: every orbit is one point
}

// ---- The molecular point-group hooks (OOD-SOLID V-FOLD-NEUTRAL step 1; user ruling 2026-11-23, D1-D3 of lane H's design note) ----------
// A molecule's ops are CARTESIAN with tau = 0 (the SymOp convention for a point group, Symmetry/SymOp.C): they act about the molecule's own
// CENTROID, the unweighted mean of its sites, which the MAKER supplies -- so a source resolved at one geometry never places the ops of
// another (V-MESH-FACTORY's rule: the source gives ops, the maker owns the geometry).  The group handled is of D2H TYPE: every op an
// orthogonal involution diagonal in one common frame (the abelian subgroup the SALC blocking uses, Symmetry::Molecule::BuildAbelianGroup).
// In that frame the one shared Lebedev set is invariant under every such op (a Lebedev rule is octahedrally invariant in its own frame),
// so the mesh is invariant with NO new rule and no site adaptation; a group beyond D2h type (a C3, a C6) is refused by name -- the full
// point group is the follow-up the design note lists, through the site-adapted builder the solid uses.

//! The molecule's ops as they act on Cartesian points: {W | c - W c}, about the centroid \a c.
static std::vector<Symmetry::SymOp> AboutCentre(const std::vector<Symmetry::SymOp>& ops, const rvec3_t& c)
{
    std::vector<Symmetry::SymOp> out;
    out.reserve(ops.size());
    for (const Symmetry::SymOp& op : ops) out.push_back({op.W, c-op.W*c, op.sigma});
    return out;
}

//! THROWS unless every op is a point-group op this build handles: tau = 0, no spin action, W orthogonal.
static void RequirePointGroupOps(const std::vector<Symmetry::SymOp>& ops, const char* who)
{
    for (size_t k=0; k<ops.size(); k++)
    {
        const Symmetry::SymOp& op=ops[k];
        if (norm(op.tau)>0.0)
            throw std::invalid_argument(std::string(who)+"(Molecule): op "+std::to_string(k)+" carries a translation; a point-group op is "
                                        "{W | 0}, acting about the molecule's centroid (Symmetry/SymOp.C)");
        if (op.sigma!=Symmetry::SpinAction::None)
            throw std::invalid_argument(std::string(who)+"(Molecule): op "+std::to_string(k)+" carries a spin action; a molecule's magnetic "
                                        "(Shubnikov) imposition is not built");
        const Matrix3D<double> P=Transpose(op.W)*op.W;
        for (size_t i=1; i<=3; i++)
            for (size_t j=1; j<=3; j++)
                if (std::abs(P(i,j)-(i==j ? 1.0 : 0.0))>1e-10)
                    throw std::invalid_argument(std::string(who)+"(Molecule): op "+std::to_string(k)+" is not orthogonal");
    }
}

//! \brief THE FRAME of a D2h-type op set, from the ops alone (no geometry).  Each op W != E is a symmetric involution, so (I+W)/2 and
//! (I-W)/2 project onto its +1 and -1 eigenspaces; a rank-1 projector names an AXIS (a C2 axis, a mirror normal).  The axes of a D2h-type
//! set are mutually orthogonal; completed to an orthonormal frame (the lab axes taken where they are free), they give the frame every
//! op is diagonal (+-1) in.  The frame may be improper: the Lebedev set is invariant under the frame's reflections too.  The frame keeps the LAB axes wherever the ops allow (each frame axis is matched to the lab
//! axis it lies closest to, sign positive), so a deck that already aligns its molecule with x, y, z gets the identity: the same mesh as
//! the plain build.  Sets \a why (and returns the identity) when an op is not diagonal in the frame: a C3, C4 or C6, outside D2h type.
static Matrix3D<double> D2hFrame(const std::vector<Symmetry::SymOp>& ops, std::string& why)
{
    std::vector<rvec3_t> axes;
    auto addAxis=[&](const rvec3_t& a)
    {
        for (const rvec3_t& b : axes) if (std::abs(a*b)>1.0-1e-8) return;   // the same axis (either sign)
        axes.push_back(a);
    };
    const Matrix3D<double> I(1,0,0, 0,1,0, 0,0,1);
    for (const Symmetry::SymOp& op : ops)
        for (double s : {+1.0, -1.0})
        {
            Matrix3D<double> Pr=I+op.W*s;   // 2x the projector onto the s eigenspace
            const double tr=Pr(1,1)+Pr(2,2)+Pr(3,3);
            if (std::abs(tr-2.0)>1e-8) continue;   // not rank 1 (an involution's projector has trace = rank, here doubled)
            rvec3_t best(0,0,0);
            for (size_t j=1; j<=3; j++)
            {
                const rvec3_t col(Pr(1,j), Pr(2,j), Pr(3,j));
                if (norm(col)>norm(best)) best=col;
            }
            if (norm(best)>1e-8) addAxis(best/norm(best));
        }
    for (size_t i=0; i<axes.size(); i++)
        for (size_t j=i+1; j<axes.size(); j++)
            if (std::abs(axes[i]*axes[j])>1e-8)
            {
                why="two of its symmetry axes are not perpendicular (a C3 or C6 family, outside D2h type)";
                return I;
            }
    if (axes.size()>3) { why="it has more than three symmetry axes (outside D2h type)"; return I; }
    // Complete with the lab axes, least aligned with the span first (Gram-Schmidt).
    const rvec3_t lab[3]={rvec3_t(1,0,0), rvec3_t(0,1,0), rvec3_t(0,0,1)};
    while (axes.size()<3)
    {
        int pick=-1; double least=1e300;
        for (int k=0; k<3; k++)
        {
            double proj=0.0; for (const rvec3_t& a : axes) proj+=(a*lab[k])*(a*lab[k]);
            if (proj<least-1e-12) { least=proj; pick=k; }
        }
        rvec3_t v=lab[pick];
        for (const rvec3_t& a : axes) v-=(a*v)*a;
        axes.push_back(v/norm(v));
    }
    // Match each lab axis to the frame axis closest to it (greedy on |cos|), sign positive: the identity when the deck is aligned.
    Matrix3D<double> F;
    std::vector<bool> used(3, false);
    for (int k=0; k<3; k++)
    {
        int best=-1; double bc=-1.0;
        for (int a=0; a<3; a++) if (!used[a] && std::abs(axes[a]*lab[k])>bc+1e-12) { bc=std::abs(axes[a]*lab[k]); best=a; }
        used[best]=true;
        const rvec3_t v=(axes[best]*lab[k]<0.0) ? -1.0*axes[best] : axes[best];
        F(1,k+1)=v.x; F(2,k+1)=v.y; F(3,k+1)=v.z;
    }
    for (size_t k=0; k<ops.size(); k++)
    {
        const Matrix3D<double> D=Transpose(F)*ops[k].W*F;
        for (size_t i=1; i<=3; i++)
            for (size_t j=1; j<=3; j++)
            {
                const double want = i!=j ? 0.0 : (D(i,i)>0.0 ? 1.0 : -1.0);
                if (std::abs(D(i,j)-want)>1e-8)
                {
                    why="op "+std::to_string(k)+" is not a reflection of the frame axes (outside D2h type: a C3, C4 or C6)";
                    return I;
                }
            }
    }
    return F;
}

//! The nearest point of \a pts within \a tol of \a r, or -1: a sorted-by-x window search (the Euclidean twin of the fold's matcher).
class EuclideanIndex
{
public:
    EuclideanIndex(const rvec3vec_t& pts, double tol) : itsPts(pts), itsTol(tol), itsByX(pts.size())
    {
        for (size_t i=0; i<pts.size(); i++) itsByX[i]=i;
        std::sort(itsByX.begin(), itsByX.end(), [&](size_t a, size_t b){return itsPts[a].x<itsPts[b].x;});
    }
    long Find(const rvec3_t& r) const
    {
        auto lo=std::lower_bound(itsByX.begin(), itsByX.end(), r.x-itsTol, [&](size_t a, double v){return itsPts[a].x<v;});
        long best=-1; double bestD2=itsTol*itsTol;
        for (auto it=lo; it!=itsByX.end() && itsPts[*it].x<=r.x+itsTol; ++it)
        {
            const rvec3_t d=itsPts[*it]-r;
            const double d2=d*d;
            if (d2<=bestD2) { best=long(*it); bestD2=d2; }
        }
        return best;
    }
private:
    const rvec3vec_t&   itsPts;
    double              itsTol;
    std::vector<size_t> itsByX;
};

//! The sites of a structure, as the polish and the geometry checks read them.
struct SitePoints
{
    std::vector<rvec3_t> R;
    std::vector<int>     Z;
    rvec3_t              c=rvec3_t(0,0,0);   //!< the unweighted centroid: the point the ops act about
};
static SitePoints SitesOf(const Structure& s)
{
    SitePoints p;
    s.ForEachSite([&p](int Z, const rvec3_t& r, bool){ p.R.push_back(r); p.Z.push_back(Z); });
    for (const rvec3_t& r : p.R) p.c+=r;
    if (!p.R.empty()) p.c=p.c/double(p.R.size());
    return p;
}

//! The nearest site of species \a Z to \a r (no tolerance: the caller judges the distance).
static size_t NearestSite(const SitePoints& p, int Z, const rvec3_t& r)
{
    size_t best=0; double bd=1e300;
    for (size_t j=0; j<p.R.size(); j++)
        if (p.Z[j]==Z && norm(p.R[j]-r)<bd) { bd=norm(p.R[j]-r); best=j; }
    return best;
}

//! \brief The D2h-type ops POLISHED against the geometry: the same group, exact to roundoff for a geometry that carries it to roundoff.
//! The ops a source hands over are DETECTED: built from axes found to the detection tolerance, so on a deck turned off the lab axes they
//! are exact only to that axis search, and a mesh carried by them would inherit its error as an orientation dependence.  So the maker
//! re-derives each frame axis from the ATOMS (the geometry it owns, V-MESH-FACTORY's rule), through an
//! op that has that axis as its one-dimensional eigenspace (eigenvalue s: a C2 axis, s = +1; a mirror normal, s = -1).  With x = R_i - c
//! and y = W x its image site's offset, x + s y = 2 (x.v) v exactly: the axis is the largest such vector.  When every such vector
//! vanishes (all atoms in the plane perpendicular to v) the axis is that plane's normal, the largest cross product of two offsets; when
//! that vanishes too (a linear molecule's free axes) nothing fixes it, and the detected axis is kept.  The axes are re-orthonormalized
//! (the fixed ones first), and the ops rebuilt as W = F D F^T from the frame's sign matrices D.  Pure: ops in, ops out, no tolerance.
static std::vector<Symmetry::SymOp> PolishedOps(const SitePoints& p, const std::vector<Symmetry::SymOp>& ops, const Matrix3D<double>& F)
{
    std::vector<Matrix3D<double>> D(ops.size());
    for (size_t k=0; k<ops.size(); k++)
    {
        const Matrix3D<double> M=Transpose(F)*ops[k].W*F;
        D[k]=Matrix3D<double>(M(1,1)>0 ? 1 : -1,0,0, 0,M(2,2)>0 ? 1 : -1,0, 0,0,M(3,3)>0 ? 1 : -1);
    }
    rvec3_t axis[3];
    bool    fixed[3]={false, false, false};
    double  scale=0.0;
    for (const rvec3_t& r : p.R) scale=std::max(scale, norm(r-p.c));
    const double zero=1e-6*std::max(scale, 1.0);   // "vanishes": far below any bond, far above roundoff
    for (int j=0; j<3; j++)
    {
        const rvec3_t approx(F(1,j+1), F(2,j+1), F(3,j+1));
        axis[j]=approx;
        for (size_t k=0; k<ops.size() && !fixed[j]; k++)
        {
            const double s=D[k](j+1,j+1);
            if (D[k]((j+1)%3+1,(j+1)%3+1)!=-s || D[k]((j+2)%3+1,(j+2)%3+1)!=-s) continue;   // axis j is not this op's 1-D eigenspace
            rvec3_t best(0,0,0);
            for (size_t i=0; i<p.R.size(); i++)
            {
                const rvec3_t x=p.R[i]-p.c;
                const rvec3_t y=p.R[NearestSite(p, p.Z[i], ops[k].W*x+p.c)]-p.c;
                const rvec3_t v=x+s*y;
                if (norm(v)>norm(best)) best=v;
            }
            if (norm(best)<=zero)   // every atom in the plane perpendicular to the axis: its normal
                for (size_t a=0; a<p.R.size(); a++)
                    for (size_t b=a+1; b<p.R.size(); b++)
                    {
                        const rvec3_t v=Cross(p.R[a]-p.c, p.R[b]-p.c);
                        if (norm(v)>norm(best)) best=v;
                    }
            if (norm(best)>zero*zero) { axis[j]=(best*approx<0.0 ? -1.0 : 1.0)*best/norm(best); fixed[j]=true; }
        }
    }
    // Gram-Schmidt, the atom-fixed axes first (they are orthogonal to roundoff already; this only removes the roundoff).
    int order[3], n=0;
    for (int j=0; j<3; j++) if ( fixed[j]) order[n++]=j;
    for (int j=0; j<3; j++) if (!fixed[j]) order[n++]=j;
    for (int a=0; a<3; a++)
    {
        rvec3_t v=axis[order[a]];
        for (int b=0; b<a; b++) v-=(axis[order[b]]*v)*axis[order[b]];
        axis[order[a]]=v/norm(v);
    }
    Matrix3D<double> P;
    for (int j=0; j<3; j++) { P(1,j+1)=axis[j].x; P(2,j+1)=axis[j].y; P(3,j+1)=axis[j].z; }
    std::vector<Symmetry::SymOp> out;
    out.reserve(ops.size());
    for (size_t k=0; k<ops.size(); k++) out.push_back({P*D[k]*Transpose(P), ops[k].tau, ops[k].sigma});
    return out;
}

//! How far the sites are from carrying \a ops about their centroid (bohr): the largest distance from an op's image of a site to the
//! nearest site of its species.
static double GeometryError(const SitePoints& p, const std::vector<Symmetry::SymOp>& ops)
{
    double d=0.0;
    for (const Symmetry::SymOp& op : ops)
        for (size_t i=0; i<p.R.size(); i++)
        {
            const rvec3_t img=op.W*(p.R[i]-p.c)+p.c;
            d=std::max(d, norm(p.R[NearestSite(p, p.Z[i], img)]-img));
        }
    return d;
}

//! The frame of \a ops and the ops polished against \a s; THROWS by name outside D2h type.
static std::pair<Matrix3D<double>, std::vector<Symmetry::SymOp>> FrameAndPolishedOps(const SitePoints& p,
                                                                                     const std::vector<Symmetry::SymOp>& ops, const char* who)
{
    std::string why;
    const Matrix3D<double> F0=D2hFrame(ops, why);
    if (!why.empty())
        throw std::invalid_argument(std::string(who)+"(Molecule, "+std::to_string(ops.size())+" ops): a mesh invariant under this group is not "
                                    "built: "+why+".  The molecular mesh is built for groups of D2h type, the abelian subgroup the SALC blocking "
                                    "uses (OOD-SOLID V-FOLD-NEUTRAL; the full point group is its follow-up)");
    std::vector<Symmetry::SymOp> polished=PolishedOps(p, ops, F0);
    std::string why2;
    const Matrix3D<double> F=D2hFrame(polished, why2);   // the polished axes, matched to the lab axes as before
    assert(why2.empty());
    return {F, std::move(polished)};
}

double PointGroupGeometryError(const Structure& s, const std::vector<Symmetry::SymOp>& ops)
{
    RequirePointGroupOps(ops, "PointGroupGeometryError");
    const SitePoints p=SitesOf(s);
    std::string why;
    (void)D2hFrame(ops, why);
    if (!why.empty()) return std::numeric_limits<double>::infinity();
    return GeometryError(p, FrameAndPolishedOps(p, ops, "PointGroupGeometryError").second);
}

qcMesh::Mesh MakeMesh(const Molecule& mol, const qcMesh::MeshParams& mp, const std::vector<Symmetry::SymOp>& ops)
{
    assert(!ops.empty());
    RequirePointGroupOps(ops, "MakeMesh");
    Symmetry::Lattice_3D::RequireClosedGroup(ops, "MakeMesh(Molecule, ops)");
    if (mp.angular==qcMesh::AngularKind::None)
        throw std::invalid_argument("MakeMesh(Molecule, ops): the radial-only set exists only on an atom (doc/Pins.md pin 30)");

    const SitePoints sites=SitesOf(mol);
    MolecularOrbitBuild sym;
    auto [frame, polished]=FrameAndPolishedOps(sites, ops, "MakeMesh");
    sym.frame=frame;

    if (mp.angRot!=0.0)
        throw std::invalid_argument("MakeMesh(Molecule, ops): the angular rule is stated ROTATED (angRot "+std::to_string(mp.angRot)+" rad), an "
                                    "orientation of its own that the group's frame would undo -- an invariant mesh is built from an unrotated "
                                    "rule (the molecule's setup keeps a rotated one in the lab frame, policy row MeshSymmetry)");
    // The angular set, oriented, must carry every op (as the frame's sign matrix).  Lebedev always does; a GaussLegendre product set
    // carries a frame reflection only when its azimuthal nodes are symmetric -- REFUSED by name, before any point is built.
    {
        const qcMesh::AngularMesh a=qcMesh::MakeAngular(mp);
        EuclideanIndex index(a.Dirs(), 1e-10);
        for (size_t k=0; k<ops.size(); k++)
        {
            const Matrix3D<double> D=Transpose(sym.frame)*ops[k].W*sym.frame;
            for (size_t i=0; i<a.size(); i++)
            {
                const long j=index.Find(D*a.Dirs()[i]);
                if (j<0 || std::abs(a.W()[size_t(j)]-a.W()[i])>1e-12*std::abs(a.W()[i]))
                    throw std::invalid_argument(std::string("MakeMesh(Molecule, ops): the ")+AngularName(mp.angular)+" angular rule of degree "
                                                +std::to_string(mp.angularDegree)+" ("+std::to_string(a.size())+" directions) does not carry the "
                                                "molecule's point group (op "+std::to_string(k)+": a reflection of the group's frame) -- To fix "
                                                "this, please state a Lebedev rule of degree 9 or more, each of which carries every op of D2h type "
                                                "(the 24- and 30-direction rules, degrees 7 and 8, do not)");
            }
        }
    }

    // The ATOM ORBITS and the edge op onto each member.  The geometry must carry the ops to the mesh match tolerance (a tenth of it, so the
    // carried points stay within it): a deck symmetric only to the detection tolerance (1e-4) cannot have an invariant mesh.
    const size_t natom=mol.GetNumAtoms();
    const std::vector<Symmetry::SymOp> opsC=AboutCentre(polished, sites.c);
    const std::vector<rvec3_t>& Rv=sites.R;
    const std::vector<int>&     Z =sites.Z;
    rvec3vec_t R(natom);
    for (size_t i=0; i<natom; i++) R[i]=Rv[i];
    const double atomTol=kMolecularMeshGeometryTol;
    {
        EuclideanIndex index(R, atomTol);
        for (size_t k=0; k<opsC.size(); k++)
            for (size_t i=0; i<natom; i++)
            {
                const long j=index.Find(opsC[k].W*R[i]+opsC[k].tau);
                if (j<0 || Z[size_t(j)]!=Z[i])
                    throw std::invalid_argument("MakeMesh(Molecule, ops): op "+std::to_string(k)+" does not map atom "+std::to_string(i)
                                                +" onto an atom of its element within "+std::to_string(atomTol)+" bohr -- the geometry does "
                                                "not carry this point group to the precision an invariant mesh needs");
            }
    }
    const Symmetry::Lattice_3D::Fold af=Symmetry::Lattice_3D::FoldPoints(Rv, opsC, atomTol);
    sym.repOf.assign(natom, 0);
    sym.edge.assign(natom, Matrix3D<double>(1,0,0, 0,1,0, 0,0,1));
    for (size_t o=0; o<af.Reps(); o++)
        for (auto [mi, oi] : af.members[o])
        {
            if (oi<0) throw std::logic_error("MakeMesh(Molecule, ops): an atom orbit member with no edge op in a closed group");
            sym.repOf[size_t(mi)]=size_t(af.repRaw[o]);
            sym.edge [size_t(mi)]=polished[size_t(oi)].W;
        }

    qcMesh::Mesh m=AtomCentredMesh(mol, mp, &sym);
    RequireSiteBlocks(m, "MakeMesh(Molecule, ops)");

    // VERIFIED, not assumed: every op maps the mesh onto itself, weight for weight.  The members are carried copies (bit-equal); what this
    // checks is each representative's own stabilizer (an atom on a mirror, at the centre), whose image points carry INDEPENDENTLY computed
    // partition weights.  Those agree to ROUNDOFF, and roundoff here grows with the distance: the Becke cell variable is a difference of two
    // distances to the point over the atoms' separation, so at distance r from the nuclei it carries ~eps r / R_min, and the cell polynomial
    // amplifies that.  Measured: ~3e-14 relative near the nuclei (MolecularPointGroupMesh tests), 9e-9 on water's outermost radial shell
    // at r = 1e6 bohr (M_Calculation.WaterLDA_eqRotated, 2026-11-23), i.e. ~1.7e-14 per unit of r / R_min.  So the spread is NORMALIZED
    // by max(1, r / R_min) before it is judged: a wrong frame or a mismatched edge op gives O(1) spreads at every r; past 1e-10 normalized
    // it is a broken invariant, not roundoff.  The normalized spread is reported.
    // The spread is that of the PARTITION, |ΔP| = |Δw| / (w_q), not |Δw| / w (lane H 2026-11-24): P ∈ [0,1] is what carries absolute
    // roundoff, and a point deep in another atom's cell has P down to 1e-26, whose RELATIVE roundoff is not bounded by anything.  Measured on
    // SF₆'s pilot-raised meshes (to 726k points): |Δw|/w reached 1.2e-10 normalized at P ~ 1e-15 and threw on a correct mesh; at P ~ 1 it
    // stays at the 1e-16..1e-14 above.  The template weight w_q = w / P, with P recomputed by the build's own qcMesh::BeckePartitionWeight.
    double rMin=std::numeric_limits<double>::infinity();
    for (size_t a=0; a<natom; a++)
        for (size_t b=a+1; b<natom; b++) rMin=std::min(rMin, norm(Rv[a]-Rv[b]));
    if (!std::isfinite(rMin)) rMin=1.0;
    double spread=0.0;
    {
        EuclideanIndex index(m.Points(), kMeshMatchTol);
        const rvec_t& w=m.Weights();
        std::vector<size_t> siteOf(m.size(), 0);
        for (size_t a=0; a<m.NSites(); a++) for (size_t i=m.SiteBegin(a); i<m.SiteEnd(a); i++) siteOf[i]=a;
        const int beckeOrder=mp.beckeOrder;
        for (size_t k=0; k<opsC.size(); k++)
            for (size_t i=0; i<m.size(); i++)
            {
                const long j=index.Find(opsC[k].W*m.Points()[i]+opsC[k].tau);
                if (j<0)
                    throw std::logic_error("MakeMesh(Molecule, ops): the built mesh is not invariant under op "+std::to_string(k)
                                           +": point "+std::to_string(i)+" has no image point");
                if (w[i]==0.0) continue;
                const double P=qcMesh::BeckePartitionWeight(Rv, siteOf[i], m.Points()[i], beckeOrder);
                if (P==0.0) continue;
                const double reach=std::max(1.0, norm(m.Points()[i]-sites.c)/rMin);
                spread=std::max(spread, std::abs(w[size_t(j)]-w[i])*P/std::abs(w[i])/reach);   // |ΔP| = |Δw| / w_q, w_q = w / P
            }
    }
    if (spread>1e-10)
    {
        std::ostringstream o; o << std::scientific << spread;
        throw std::logic_error("MakeMesh(Molecule, ops): the weights differ across an orbit by "+o.str()+" relative (per unit of r / R_min) "
                               "-- not roundoff");
    }
    qchem::report::Set("orbitWeightSpread", spread);
    {
        std::ostringstream fs;
        fs << std::setprecision(6);
        for (size_t i=1; i<=3; i++) fs << (i>1 ? "; " : "") << sym.frame(i,1) << ' ' << sym.frame(i,2) << ' ' << sym.frame(i,3);
        qchem::report::Set("pointGroupOps", (long)ops.size());
        qchem::report::Set("frame", fs.str());   // rows of the matrix whose columns are the frame axes
    }
    return m;
}

qcMesh::FoldedMesh FoldOrbits(const Molecule& mol, qcMesh::Mesh mesh, const std::vector<Symmetry::SymOp>& ops)
{
    RequirePointGroupOps(ops, "FoldOrbits");
    const SitePoints sites=SitesOf(mol);   // the ops the mesh was built under: polished against the same geometry
    const std::vector<Symmetry::SymOp> opsC=AboutCentre(FrameAndPolishedOps(sites, ops, "FoldOrbits").second, sites.c);
    auto fold=[&]{ qchem::report::Timed timed("setup: XC mesh orbit fold (FoldMesh)");
                   std::vector<rvec3_t> pts(mesh.size());
                   for (size_t i=0; i<mesh.size(); i++) pts[i]=mesh.Points()[i];
                   return Symmetry::Lattice_3D::FoldPoints(pts, opsC, kMeshMatchTol); }();
    // The SITE ACTION: which atom each op carries each atom onto.  The site blocks are the atoms in ForEachSite order (RequireSiteBlocks in
    // MakeMesh), and MakeMesh already refused a geometry an op does not map onto itself, so every image is found.
    qcMesh::FoldedMesh::SiteAction act;
    {
        const size_t natom=sites.R.size();
        if (mesh.NSites()!=natom)
            throw std::logic_error("FoldOrbits(Molecule): the mesh has "+std::to_string(mesh.NSites())+" site blocks for "+std::to_string(natom)
                                   +" atoms -- not a mesh MakeMesh(Molecule, ops) built");
        rvec3vec_t R(natom);
        for (size_t i=0; i<natom; i++) R[i]=sites.R[i];
        EuclideanIndex index(R, kMolecularMeshGeometryTol);
        for (const Symmetry::SymOp& op : opsC)
        {
            std::vector<size_t> img(natom);
            for (size_t i=0; i<natom; i++)
            {
                const long j=index.Find(op.W*R[i]+op.tau);
                if (j<0 || sites.Z[size_t(j)]!=sites.Z[i])
                    throw std::logic_error("FoldOrbits(Molecule): an op does not map atom "+std::to_string(i)+" onto an atom of its element, "
                                           "though MakeMesh accepted the geometry");
                img[i]=size_t(j);
            }
            act.W.push_back(op.W);
            act.image.push_back(std::move(img));
        }
    }
    std::cout << "[Becke grid] molecular point group (" << ops.size() << " ops): invariant mesh " << mesh.size() << " points in "
              << fold.Reps() << " orbits; rho star-averaged each iteration" << std::endl;
    return {std::make_shared<const qcMesh::Mesh>(std::move(mesh)), std::move(fold), std::move(act)};
}

} // namespace qchem
