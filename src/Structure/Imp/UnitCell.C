// File: Structure/UnitCell.C  Unit cell for a lattice.  Symbol/units conventions: see Lattice.C.
module;
#include <iomanip>
#include <cassert>
#include <cstdlib>   // std::getenv (the GPW_OMP_THREADS cap on the Becke mesh build)
#include <iostream>
#include <string>
#include <vector>

module qchem.UnitCell;
import qchem.Math;
import qchem.Structure;      // Atom (AddAtom inserts atoms given in fractional coordinates)
import qchem.Vector3D;       // norm(rvec3_t)
import qchem.Mesh.Product;   // ProductMesh, MakeRadial, MakeAngular (the single-centre template)
import qchem.Mesh.Builder;   // MeshBuilder (efficient incremental accumulation)
import qchem.SymmetrizeMesh; // MakeInvariantAngularMesh + FoldPointsPeriodic (the §6a W2b site-adapted mesh)
import qchem.Reporting;      // grids.becke report addendum (EmitAt; inert when no run is open)

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

// One competitor image atom R_b + L for the Becke partition at one mesh point.
struct BeckeImage
{
    rvec3_t R;        // image position
    double  d;        // |r - R|
    bool    inPSet;   // cell-function bound >= eps: P is computed and enters the denominator
};

//! \a angPerAtom (optional) = per-atom angular sets -- the §6a W2b SITE-ADAPTED path (rep atoms carry
//! site-group-invariant sets, partners the op-rotated copies; the ops-taking CreateIntegrationMesh
//! overload builds them).
//! Null = the historical shared single-orientation angular grid.
qcMesh::Mesh MakePeriodicBeckeMesh(const UnitCell& cell, const qcMesh::MeshParams& mp,
                                   const std::vector<qcMesh::AngularMesh>* angPerAtom = nullptr)
{
    qchem::report::Timed timed("setup: becke mesh build");
    using qcMesh::BeckeCutoff;
    const int    k  =mp.beckeOrder;
    // Magnitude screen: the weights are eps-converged.  1e-8 (not tighter) on PURPOSE: the tail bounds are
    // worst-case-COLLINEAR (attained only on the exact line through two images), so the true neglected mass
    // is far below eps -- while the quadrature error itself sits at ~1e-4.  Tightening to 1e-10 was measured
    // to inflate the competitor sets (and the O(|P-set| x |images|) partition cost) several-fold for sparse
    // cells with no accuracy gain at the quadrature level.
    const double eps=1e-8;
    assert(k>=0);

    std::vector<rvec3_t> R;
    for (auto a : cell) R.push_back(a->itsR);
    const size_t natom=R.size();
    assert(natom>0);

    // Minimum interplanar spacing d_i = 1/sqrt((M^-1)_ii): a point in cell n0 is at least s*dPlane from
    // anything at Chebyshev cell distance > s (s whole cell layers in between) -- the shell-tail floor.
    const Matrix3D<double>& A=cell.GetCellMatrix();
    Matrix3D<double> Minv=Invert(Transpose(A)*A);
    const double dPlane=1.0/sqrt(max(Minv(1,1),max(Minv(2,2),Minv(3,3))));

    qcMesh::RadialMesh  rad=qcMesh::MakeRadial(mp);   // one single-centre template reused (ShiftOrigin per atom)
    qcMesh::AngularMesh ang=qcMesh::MakeAngular(mp);

    qcMesh::MeshBuilder out;
    for (size_t ia=0; ia<natom; ia++)
    {
        qcMesh::Mesh am=qcMesh::ProductMesh(rad, angPerAtom ? (*angPerAtom)[ia] : ang);
        am.ShiftOrigin(R[ia]);
        const rvec3vec_t& pts=am.Points();
        const rvec_t&     wts=am.Weights();
        // The per-point Becke partitions are INDEPENDENT (~1 ms each -- the eps-converged competitor
        // image series is the whole build cost, measured linear at any angular count), and each result
        // lands in its own indexed slot consumed in index order below -- so the threaded build is
        // BIT-IDENTICAL to the serial one at any thread count (no cross-point reductions, unlike the
        // GPW pair loops).  Hence parallel by DEFAULT under QCHEM_OPENMP; GPW_OMP_THREADS, when set,
        // is honoured as the thread CAP (user request: ONE knob governs all the heavy loops) --
        // purely resource control, never a numerics knob.
        const size_t nq=am.size();
        std::vector<char> keep(nq, 0);
        rvec3vec_t        kpt(nq);
        rvec_t            kwt(nq);
        auto partition=[&](size_t q)
        {
            std::vector<BeckeImage> im;   // competitor images for this point (thread-local)
            const rvec3_t& r=pts[q];
            // The point's own cell n0 (shells are centred on it, so the s*dPlane floor holds).
            rvec3_t f=cell.ToFractional(r);
            const int n0x=int(floor(f.x)), n0y=int(floor(f.y)), n0z=int(floor(f.z));

            long   target=-1;                 // index of THIS grid's own centre (atom ia, L=0) in im
            double dNear=1e300;
            auto AddShell=[&](int s)
            {
                for (int i=-s; i<=s; i++)
                    for (int j=-s; j<=s; j++)
                        for (int l=-s; l<=s; l++)
                        {
                            if (max(abs(i),max(abs(j),abs(l)))!=s) continue;   // surface of the cube only
                            const int nx=n0x+i, ny=n0y+j, nz=n0z+l;
                            rvec3_t L=cell.ToCartesian(rvec3_t(nx,ny,nz));
                            for (size_t b=0; b<natom; b++)
                            {
                                rvec3_t Rb=R[b]+L;
                                double  d=norm(r-Rb);
                                if (b==ia && nx==0 && ny==0 && nz==0) target=long(im.size());
                                im.push_back({Rb,d,false});
                                if (d<dNear) dNear=d;
                            }
                        }
            };
            AddShell(0);
            AddShell(1);
            // Early drop of far radial-tail points: dNear over shells 0..1 only DECREASES as more shells
            // arrive, and s((da-dn)/(da+dn)) is increasing in dn, so this is a valid upper bound on the
            // point's own cell function.
            const double da=norm(r-R[ia]);
            if (BeckeCutoff((da-dNear)/(da+dNear),k)<eps) return;

            for (int s=2; ; s++)
            {
                const double Dlow=(s-1)*dPlane;   // distance floor for anything at cell distance >= s
                const double newPB=BeckeCutoff((Dlow-dNear)/(Dlow+dNear),k);
                if (newPB<eps)
                {
                    double facB=0;
                    for (auto& b : im)
                    {
                        double pb=BeckeCutoff((b.d-dNear)/(b.d+dNear),k);
                        b.inPSet=(pb>=eps);
                        if (!b.inPSet) continue;
                        double dev=pb*(1.0-BeckeCutoff((b.d-Dlow)/(b.d+Dlow),k));
                        if (dev>facB) facB=dev;
                    }
                    if (facB<eps) break;          // both tails below eps: the series is converged
                }
                AddShell(s);
                assert(s<100);                    // physically unreachable; guards a runaway series
            }

            // If the point's own centre was never gathered, its cell function is below the unseen-image
            // bound (< eps) -- the partition assigns it ~0 weight.
            if (target<0 || !im[target].inPSet) return;

            double sum=0, Pt=0;
            for (size_t i=0; i<im.size(); i++)
            {
                if (!im[i].inPSet) continue;
                double P=1.0;
                for (size_t j=0; j<im.size() && P>0; j++)
                {
                    if (j==i) continue;
                    double Rij=norm(im[i].R-im[j].R);
                    double mu = (Rij>0.0) ? (im[i].d-im[j].d)/Rij : 0.0;  // coincident atoms -> mu=0
                    P*=BeckeCutoff(mu,k);
                    if (P<1e-300) P=0.0;          // underflow guard, not a screen
                }
                sum+=P;
                if (long(i)==target) Pt=P;
            }
            assert(sum>0);
            double w=Pt/sum;
            if (w>0)   // wrapped into the home cell; slot-indexed so the threaded build stays ordered
            {
                keep[q]=1;
                kpt[q]=r-cell.ToCartesian(rvec3_t(n0x,n0y,n0z));
                kwt[q]=wts[q]*w;
            }
        };
#ifdef QCHEM_OPENMP
        static const int cap=[]{ const char* s=std::getenv("GPW_OMP_THREADS"); return s?std::atoi(s):0; }();
        if (cap>0)
        {
            #pragma omp parallel for schedule(dynamic, 8) num_threads(cap)
            for (size_t q=0; q<nq; q++) partition(q);
        }
        else
        {
            #pragma omp parallel for schedule(dynamic, 8)
            for (size_t q=0; q<nq; q++) partition(q);
        }
#else
        for (size_t q=0; q<nq; q++) partition(q);
#endif
        for (size_t q=0; q<nq; q++) if (keep[q]) out.Append(kpt[q], kwt[q]);
        // GPW_BECKE_ATOMS: the PER-ATOM breakdown.  Sum(w) is the partition's share of unity assigned to
        // this atom -- for two crystallographically EQUIVALENT atoms it must be identical, and so must the
        // kept/dropped counts.  The aggregate "[Becke grid] ... dropped N tail pts" line cannot show that,
        // which is why a site-dependent mesh defect could hide behind it (MnO 2026-08-11: the moment dies
        // on the atom at the cell CORNER and survives on the one at the cell CENTRE, and a rigid
        // translation of the whole crystal -- an exact symmetry -- moves the answer by 56 Ha).
        if (std::getenv("GPW_BECKE_ATOMS"))
        {
            size_t nk=0; double wsum=0.0;
            for (size_t q=0; q<nq; q++) if (keep[q]) { ++nk; wsum+=kwt[q]; }
            std::cout<<"[Becke atom] ia="<<ia<<" R=("<<R[ia].x<<","<<R[ia].y<<","<<R[ia].z<<")"
                     <<" pts="<<nq<<" kept="<<nk<<" dropped="<<(nq-nk)
                     <<" Sum(w)="<<std::setprecision(10)<<wsum<<std::endl;
        }
    }
    qcMesh::Mesh mesh=out.take();

    // Announce the grid -- the console line (the [GPW grid] family: one line, run start) is HOW a user
    // tells the Becke XC quadrature is in play at all; the report addendum lands in the `grids` section
    // next to the uniform ladder whenever a report run is open (inert otherwise, e.g. bare mesh tests).
    auto AngularName=[](qcMesh::AngularKind a)
    { switch (a) {
        case qcMesh::AngularKind::Lebedev:         return "Lebedev";
        case qcMesh::AngularKind::GaussLegendre: return "GaussLegendre"; }
      return "?"; };
    auto RadialName=[](qcMesh::RadialKind r)
    { switch (r) {
        case qcMesh::RadialKind::MHL:    return "MHL";
        case qcMesh::RadialKind::Log:    return "Log";
        case qcMesh::RadialKind::Linear: return "Linear"; }
      return "?"; };
    size_t nDirs=0;
    if (angPerAtom) for (const auto& a : *angPerAtom) nDirs=max(nDirs, a.size());
    else            nDirs=ang.size();
    size_t total=0;
    if (angPerAtom) for (const auto& a : *angPerAtom) total+=size_t(rad.size())*a.size();
    else            total=natom*size_t(rad.size())*ang.size();
    const size_t kept=mesh.size();
    std::cout<<"[Becke grid] natom="<<natom<<" radial="<<RadialName(mp.radial)<<" nR="<<mp.nRadial
             <<" alpha="<<mp.mhl_alpha<<" angular="<<(angPerAtom?"SITE-ADAPTED":AngularName(mp.angular))
             <<" L="<<mp.angularDegree<<" ("<<nDirs<<" dirs"<<(angPerAtom?" max/atom":"")<<") beckeOrder="<<k
             <<" points="<<kept<<" (dropped "<<total-kept<<" tail pts)"<<std::endl;
    qchem::report::EmitAt("grids", "becke", {
        {"natom", (long)natom}, {"radial",RadialName(mp.radial)}, {"nRadial", mp.nRadial},
        {"mhl_alpha", mp.mhl_alpha},
        {"angular", angPerAtom ? "SiteAdapted" : AngularName(mp.angular)}, {"L", mp.angularDegree},
        {"nDirs", (long)nDirs},
        {"beckeOrder", k}, {"points", (long)kept}, {"dropped", (long)(total-kept)}, {"eps", eps}});
    return mesh;
}

} //anon

// §6a W2b: the SITE-ADAPTED periodic Becke mesh.  Per atom ORBIT of the imposed ops, the representative
// atom carries a site-group-INVARIANT angular set (MakeInvariantAngularMesh under its CARTESIAN stabilizer
// A·W·A⁻¹) and every symmetry partner carries the op-rotated copy (the atom-fold's edge op) -- so the whole
// mesh maps onto itself under every op BY CONSTRUCTION (the T2 precondition), with the fuzzy-Voronoi
// partition computed on the final point set (a geometric function of the atom distances, so the weights
// inherit the invariance).
qcMesh::Mesh UnitCell::CreateIntegrationMesh(const qcMesh::MeshParams& mp,
                                             const std::vector<Symmetry::SymOp>& ops) const
{
    namespace SL = Symmetry::Lattice_3D;
    assert(!ops.empty());

    // UNIFORM kind: there is no site-adapted construction to do -- build the generic mesh and group-average
    // it onto the op orbits (W1; a no-op dedup when the grid is already commensurate).  This branch used to
    // sit in the CALLER (GPW_IBS's XC-quadrature factory), which meant a BASIS was switching on cellKind to
    // decide how a STRUCTURE builds its own mesh.  The structure owns the mesh; it owns the strategy.
    if (mp.cellKind!=qcMesh::UnitCellKind::Becke)
        return MakeInvariant(CreateIntegrationMesh(mp), GetCellMatrix(), ops, kMeshMatchTol);

    // BECKE kind: site-adapted, invariant by construction (W2b) -- the rest of this function.

    std::vector<rvec3_t> F;                         // fractional atom sites
    for (auto a : *this) F.push_back(ToFractional(a->itsR));
    SL::Fold af = SL::FoldPointsPeriodic(F, ops, 1e-6);   // atom orbits + the rep->partner edge ops

    const Matrix3D<double>& A = GetCellMatrix();
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
        const rvec3_t r0 = ToCartesian(fr);
        std::vector<std::pair<double,rvec3_t>> nb;
        double dmin = 1e300;
        for (auto a : *this)
            for (int i = -1; i <= 1; ++i)
                for (int j = -1; j <= 1; ++j)
                    for (int k = -1; k <= 1; ++k)
                    {
                        rvec3_t dR = a->itsR + ToCartesian(rvec3_t(i,j,k)) - r0;
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
    for (size_t o = 0; o < af.Reps(); ++o)
    {
        const rvec3_t& fr = F[af.repRaw[o]];
        std::vector<Matrix3D<double>> siteOps;
        for (const auto& op : ops) if (fixes(op, fr)) siteOps.push_back(cart(op.W));
        qcMesh::Mesh aq = MakeInvariantAngularMesh(siteOps, mp.angularDegree, bondDirs(fr));
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
    qcMesh::Mesh mesh = MakePeriodicBeckeMesh(*this, mp, &ang);

    // The raw point set is op-covariant BY CONSTRUCTION, but the builder's eps-tail DROP decisions
    // (the <eps screens and the w>0 keep) are computed from bit-DIFFERENT rotated distances, so a
    // borderline far-tail point can be kept on one atom while its orbit partner is dropped -- measured
    // at the production recipe (L=29, ~70k raw points).  Cure: drop every orbit-INCOMPLETE point
    // (complete <=> |orbit| x |site stabilizer| == |ops|).  Only eps-borderline points -- weight ~eps
    // or below -- can be incomplete, so this stays inside the eps-converged-series contract.
    SL::Fold f = SL::FoldPointsPeriodic([&]{ std::vector<rvec3_t> fp;
        for (size_t i = 0; i < mesh.size(); ++i) fp.push_back(ToFractional(mesh.Points()[i]));
        return fp; }(), ops, 1e-8);
    qcMesh::MeshBuilder out;
    size_t nDropped = 0;
    for (size_t o = 0; o < f.Reps(); ++o)
    {
        rvec3_t fr = ToFractional(mesh.Points()[f.repRaw[o]]);
        size_t  stab = 0;
        for (const auto& op : ops) if (fixes(op, fr)) ++stab;
        if (f.members[o].size()*stab == ops.size())
            for (auto [mi, oi] : f.members[o]) out.Append(mesh.Points()[mi], mesh.Weights()[mi]);
        else nDropped += f.members[o].size();
    }
    if (nDropped>0)
        std::cout<<"[Becke grid] orbit-consistency: dropped "<<nDropped
                 <<" eps-borderline points whose orbit partners were tail-dropped"<<std::endl;
    return out.take();
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
qcMesh::Mesh UnitCell::CreateIntegrationMesh(const qcMesh::MeshParams& mp) const
{
    if (mp.cellKind==qcMesh::UnitCellKind::Becke) return MakePeriodicBeckeMesh(*this,mp);

    // Points per axis: DERIVE from a physical grid cutoff when mp.eCut>0 (the GPW / Nyquist path), else the
    // manual mp.nUniform.  The Nyquist relation itself is qcMesh::UniformDivisions -- named once beside the
    // two MeshParams fields it relates, because the resolution CHECKS (is this mesh fine enough for the
    // basis? which grid is cheaper?) need the same mapping and must not re-derive it.
    const int n = mp.eCut>0.0 ? qcMesh::UniformDivisions(GetMaximumCellEdge(), mp.eCut) : mp.nUniform;
    assert(n>0);
    const double w=GetCellVolume()/(double(n)*n*n);   // equal weight Omega/n^3
    rvec3vec_t R(size_t(n)*n*n);
    rvec_t     W(size_t(n)*n*n);
    size_t q=0;
    for (int i=0; i<n; i++)
        for (int j=0; j<n; j++)
            for (int k=0; k<n; k++,q++)
            {
                R[q]=ToCartesian(rvec3_t((i+0.5)/n, (j+0.5)/n, (k+0.5)/n));   // r = A f, fractional midpoints
                W[q]=w;
            }

    // Announce, like the Becke branch -- so a run's console/report always says WHICH cell quadrature is in
    // play.  DEDUPED on the (params, geometry) identity: several PP terms build this same mesh in one run
    // (PP_Local, PP_NonLocal, the GPW PP quadrature) and identical repeats would read as a bug.  Setup is
    // single-threaded (see qchem.Reporting), so the static is safe.
    static std::string lastID;
    std::string id=mp.ID()+"/n"+std::to_string(n)+"/O"+std::to_string(GetCellVolume());
    if (id!=lastID)
    {
        lastID=id;
        EmitUniformGridReport("unitCellUniform", vec3_t<int>(n,n,n), mp.eCut);
    }
    return qcMesh::Mesh(std::move(R), std::move(W));
}

// The ONE uniform-cell-grid announcement (see the interface doc): console line + grids.<key> entry with
// kind, N, points, dr_i = |a_i|/N_i -- edge lengths through the metric (GetDistance on fractional unit
// vectors) -- plus eCut when the caller derived N from one.
void UnitCell::EmitUniformGridReport(const std::string& key, const vec3_t<int>& N, double eCut) const
{
    assert(N.x>0 && N.y>0 && N.z>0);
    const double dx=GetDistance(rvec3_t(1,0,0))/N.x,
                 dy=GetDistance(rvec3_t(0,1,0))/N.y,
                 dz=GetDistance(rvec3_t(0,0,1))/N.z;
    const long npts=long(N.x)*N.y*N.z;
    std::cout<<"[uniform grid] "<<key<<": N=("<<N.x<<","<<N.y<<","<<N.z<<") ("<<npts<<" pts)"
             <<" dr=("<<dx<<","<<dy<<","<<dz<<") a.u.";
    if (eCut>=0.0) std::cout<<" eCut="<<eCut<<" Ha";
    std::cout<<std::endl;
    report::json j={{"kind","Uniform"}, {"N",{N.x,N.y,N.z}}, {"points",npts}, {"dr",{dx,dy,dz}}};
    if (eCut>=0.0) j["eCut"]=eCut;
    qchem::report::EmitAt("grids", key, std::move(j));
}

//  Build the cell matrix A (columns = lattice vectors a₁,a₂,a₃) from the cell
//  lengths a,b,c and angles α,β,γ (radians), in the standard orientation
//  a₁∥x, a₂ in the xy-plane.  Then M = AᵀA reproduces the usual metric tensor
//  \f$ M_{ij}=a_i\cdot a_j \f$.
static Matrix3D<double> CellMatrix(double a, double b, double c, double α, double β, double γ)
{
    double cα=cos(α), cβ=cos(β), cγ=cos(γ), sγ=sin(γ);
    double V1=sqrt(1.0 - cα*cα - cβ*cβ - cγ*cγ + 2.0*cα*cβ*cγ); //= cell volume / (abc)
    return Matrix3D<double>(
        a  , b*cγ, c*cβ,
        0.0, b*sγ, c*(cα-cβ*cγ)/sγ,
        0.0, 0.0 , c*V1/sγ);
}

UnitCell::UnitCell(const Matrix3D<double>& A)
    : itsA(A)
    , itsM(Transpose(A)*A)
{}

UnitCell::UnitCell(double a, double b, double c, double α, double β, double γ)
    : UnitCell(CellMatrix(a,b,c,Rad(α),Rad(β),Rad(γ)))
{}

UnitCell::UnitCell(double a) : UnitCell(a,a,a,90,90,90) {}

UnitCell UnitCell::MakeReciprocalCell() const
{
    //  Solid-state convention bᵢ·aⱼ = 2π δᵢⱼ, i.e. B = 2π A⁻ᵀ.  Hence the
    //  reciprocal metric is BᵀB = (2π)² M⁻¹ (the 2π lives in B, not M).
    return UnitCell(2.0*Pi*Transpose(Invert(itsA)));
}

rvec3_t UnitCell::ToCartesian(const rvec3_t& f) const
{
    return itsA*f; // r = A f
}

rvec3_t UnitCell::ToFractional(const rvec3_t& r) const
{
    return Invert(itsA)*r; // f = A^-1 r
}

void UnitCell::AddAtom(int Z, const rvec3_t& f, bool spinFlip)
{
    Atom* a = new Atom(Z, ToCartesian(f)); // store in Cartesian a.u. (r = A f)
    a->itsSpinFlip = spinFlip;             // collinear -m sublattice bit (spin-SAD assembly; see Structure.C)
    Insert(a);
}

// FCC primitive cell: columns are the half-face-diagonal lattice vectors (h = a/2); |det A| = a^3/4.
FCCUnitCell::FCCUnitCell(double a)
    : UnitCell(Matrix3D<double>(0.0, a/2, a/2,
                                a/2, 0.0, a/2,
                                a/2, a/2, 0.0))
{}

double UnitCell::GetCellVolume() const
{
    return fabs(Determinant(itsA)); // |det A|
}

double UnitCell::GetMinimumCellEdge() const
{
    double m=itsM(1,1); // Mᵢᵢ = |aᵢ|²
    if (itsM(2,2)<m) m=itsM(2,2);
    if (itsM(3,3)<m) m=itsM(3,3);
    return sqrt(m);
}

double UnitCell::GetMaximumCellEdge() const
{
    double m=itsM(1,1); // Mᵢᵢ = |aᵢ|²
    if (itsM(2,2)>m) m=itsM(2,2);
    if (itsM(3,3)>m) m=itsM(3,3);
    return sqrt(m);
}

Vector3D<int> UnitCell::GetNumCells(double MaxDistance) const
{
    assert(MaxDistance>0);
    //  Cells per axis needed to cover a sphere of radius MaxDistance.  Uses the
    //  interplanar spacing dᵢ = 1/√(M⁻¹)ᵢᵢ (NOT the edge length |aᵢ|), so that
    //  oblique cells are not under-counted.
    Matrix3D<double> Minv=Invert(itsM);
    return Vector3D<int>(
        (int)ceil(MaxDistance*sqrt(Minv(1,1))),
        (int)ceil(MaxDistance*sqrt(Minv(2,2))),
        (int)ceil(MaxDistance*sqrt(Minv(3,3))));
}

double UnitCell::GetDistance(const rvec3_t& f) const
{
    return sqrt(f*itsM*f); // ‖A f‖ for fractional f
}

std::vector<vec3_t<int>> UnitCell::CellsInSphere(double MaxDistance) const
{
    assert(MaxDistance>0);
    std::vector<vec3_t<int>> ret;
    vec3_t<int> nc=GetNumCells(MaxDistance);
    vec3_t<int> n;
    for (n.x=-nc.x; n.x<=nc.x; n.x++)
        for (n.y=-nc.y; n.y<=nc.y; n.y++)
            for (n.z=-nc.z; n.z<=nc.z; n.z++)
                if (GetDistance(n)<=MaxDistance) ret.push_back(n);
    return ret;
}

std::ostream& UnitCell::Write(std::ostream& os) const
{
    double a=sqrt(itsM(1,1)), b=sqrt(itsM(2,2)), c=sqrt(itsM(3,3));
    double α=acos(itsM(2,3)/(b*c))/Pi*180;
    double β=acos(itsM(1,3)/(a*c))/Pi*180;
    double γ=acos(itsM(1,2)/(a*b))/Pi*180;
    os << "(a,b,c)=(" << a << "," << b << "," << c << "), "
       << "(α,β,γ)=(" << α << "," << β << "," << γ << ")°  ";
    Molecule::Write(os);
    return os;
}

} // namespace qchem