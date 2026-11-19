// File: BeckePartition.C  The molecular Becke fuzzy-Voronoi partition weight, and its GEOMETRY GRADIENT.
//
// One home for the weight a molecular (or atom) mesh point carries, so the force's grid response differentiates exactly the weight the
// mesh was built with (doc/GeometryOptimizationPlan.md §4, the XC row's part (b); lane B, 2026-11-07).  qcStructure's MakeMolecularMesh
// calls BeckePartitionWeight; the gradient lives beside it, in qcMesh, because lane S's V-NONLINEAR-MESH takes the mesh away from the
// Hamiltonian term: whoever owns the mesh then, the weight's derivative is a property of the mesh.
//
// THE PARTITION (A. D. Becke, J. Chem. Phys. 88, 2547 (1988)).  For atoms at R_i, a point r owned by atom A:
//     w_A(r) = P_A(r) / Σ_i P_i(r),   P_i(r) = Π_{j≠i} s(μ_ij),   μ_ij = (|r − R_i| − |r − R_j|) / |R_i − R_j|,
// with s = BeckeCutoff (the iterated polynomial).  No atomic-size adjustment (the production partition has none).
//
// THE GRADIENT is Enzyme's (the user's ruling, 2026-11-08: derivatives are Enzyme's, not hand-coded; doc/GeometryOptimizationPlan.md, the
// AD verdict): BeckePartitionWeight is a derivative LEAF, differentiated in place in this unit, with the point riding with its owner
// (r = R_A + a fixed offset).  The XC force's grid-response term is Σ_g (w_rad w_ang)_g e_xc(ρ_g) ∂w_g/∂R_D (Johnson, Gill and Pople,
// J. Chem. Phys. 98, 5612 (1993)), and translation invariance makes Σ_D ∂w/∂R_D = 0 exactly.  Enzyme is
// required (user ruling 2026-11-09), so the gradient never refuses.  The hand derivation it replaced (2026-11-07) is the unit test's ORACLE (src/Mesh/tests/BeckePartitionTests.C), not a
// production path.
module;
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <stdexcept>
#include <vector>
#include <enzyme/enzyme>
export module qchem.Mesh.BeckePartition;
export import qchem.Mesh;      // BeckeCutoff (inline: a leaf this unit differentiates through), rvec3_t

export namespace qchem::qcMesh
{

//! \brief The partition weight w_owner(r) (the formula in the file header), at atoms \a R, Becke order \a k.  One atom: 1.  Coincident atoms
//! (R_ij = 0) give μ = 0, as the mesh build always has.
double BeckePartitionWeight(const std::vector<rvec3_t>& R, size_t owner, const rvec3_t& r, int k);

//! \brief ∂w_owner(r)/∂R_D for every atom D, with the point riding with its owner (r = R_owner + a fixed offset), by reverse-mode AD of
//! \c BeckePartitionWeight.  Sums to zero over D.
std::vector<rvec3_t> BeckePartitionGradient(const std::vector<rvec3_t>& R, size_t owner, const rvec3_t& r, int k);

// ---- THE PERIODIC PARTITION (moved here from qcStructure's MakePeriodicBeckeMesh, lane B 2026-11-19, so its BUILD and its DERIVATIVE are
// one computation).  For a cell atom ia's template point r, the competitors are every periodic image R_b + L of the cell atoms, gathered in
// expanding Chebyshev cell shells around r's own cell until both tails are provably below eps (the magnitude screen, never a radius:
// MakePeriodicBeckeMesh's header has the bounds).  The P-set (cell-function bound >= eps) enters the denominator sum; the products run over
// every gathered image.  w = P_target / Σ_{P-set} P_i, or 0 when the point is a dropped tail.

//! One competitor image atom R_b + L at one mesh point.
struct BeckeImage
{
    rvec3_t R;        //!< image position
    double  d;        //!< |r - R|
    bool    inPSet;   //!< cell-function bound >= eps: P is computed and enters the denominator
    size_t  atom;     //!< the cell atom b it images
};

//! GPW_BECKE_COUNT=1: the per-point COST CENSUS of the partition (qcStructure's MakePeriodicBeckeMesh reads it; the counters' meaning is there).
struct BeckePointCost
{
    unsigned im=0, pset=0, shells=0, retest=0, final_=0, bound=0, dead=0, deadWork=0;
};

//! \brief The periodic Becke partition weight of cell atom \a ia's template point \a r (unwrapped: \a r = R[ia] + its template offset), over
//! the images of the cell atoms \a R in \a cell (a type with \c ToFractional / \c ToCartesian), \a dPlane the cell's minimum interplanar
//! spacing.  Fills \a im (every gathered competitor, P-set marked) and \a target (this point's own centre, R[ia] at L = 0, or -1), and
//! returns w (0 for a dropped tail).  \a census: count the inner loop's work into \a pc.
template <class Cell> double PeriodicBeckeWeight(const Cell& cell, const std::vector<rvec3_t>& R, double dPlane, size_t ia, const rvec3_t& r,
                                                 int k, double eps, std::vector<BeckeImage>& im, long& target, BeckePointCost& pc, bool census)
{
    im.clear();
    target=-1;
    const size_t natom=R.size();
    // The point's own cell n0 (shells are centred on it, so the s*dPlane floor holds).
    const rvec3_t f=cell.ToFractional(r);
    const int n0x=int(std::floor(f.x)), n0y=int(std::floor(f.y)), n0z=int(std::floor(f.z));
    double dNear=1e300;
    auto AddShell=[&](int s)
    {
        for (int i=-s; i<=s; i++)
            for (int j=-s; j<=s; j++)
                for (int l=-s; l<=s; l++)
                {
                    if (std::max(std::abs(i),std::max(std::abs(j),std::abs(l)))!=s) continue;   // surface of the cube only
                    const int nx=n0x+i, ny=n0y+j, nz=n0z+l;
                    rvec3_t L=cell.ToCartesian(rvec3_t(nx,ny,nz));
                    for (size_t b=0; b<natom; b++)
                    {
                        rvec3_t Rb=R[b]+L;
                        double  d=norm(r-Rb);
                        if (b==ia && nx==0 && ny==0 && nz==0) target=long(im.size());
                        im.push_back({Rb,d,false,b});
                        if (d<dNear) dNear=d;
                    }
                }
    };
    AddShell(0);
    AddShell(1);
    // Early drop of far radial-tail points: dNear over shells 0..1 only DECREASES as more shells arrive, and s((da-dn)/(da+dn)) is
    // increasing in dn, so this is a valid upper bound on the point's own cell function.
    const double da=norm(r-R[ia]);
    pc.bound++;
    if (BeckeCutoff((da-dNear)/(da+dNear),k)<eps) return 0.0;

    for (int s=2; ; s++)
    {
        const double Dlow=(s-1)*dPlane;   // distance floor for anything at cell distance >= s
        pc.bound++; pc.shells=unsigned(s);
        const double newPB=BeckeCutoff((Dlow-dNear)/(Dlow+dNear),k);
        if (newPB<eps)
        {
            double facB=0;
            for (auto& b : im)
            {
                double pb=BeckeCutoff((b.d-dNear)/(b.d+dNear),k);
                pc.retest++;
                b.inPSet=(pb>=eps);
                if (!b.inPSet) continue;
                double dev=pb*(1.0-BeckeCutoff((b.d-Dlow)/(b.d+Dlow),k));
                pc.retest++;
                if (dev>facB) facB=dev;
            }
            if (facB<eps) break;          // both tails below eps: the series is converged
        }
        AddShell(s);
        assert(s<100);                    // physically unreachable; guards a runaway series
    }

    // If the point's own centre was never gathered, its cell function is below the unseen-image bound (< eps) -- the partition assigns it
    // ~0 weight.
    if (target<0 || !im[target].inPSet) return 0.0;

    double sum=0, Pt=0;
    for (size_t i=0; i<im.size(); i++)
    {
        if (!im[i].inPSet) continue;
        pc.pset++;
        double P=1.0;
        unsigned jThis=0;
        for (size_t j=0; j<im.size() && P>0; j++)
        {
            if (j==i) continue;
            double Rij=norm(im[i].R-im[j].R);
            double mu = (Rij>0.0) ? (im[i].d-im[j].d)/Rij : 0.0;  // coincident atoms -> mu=0
            P*=BeckeCutoff(mu,k);
            if (census) { pc.final_++; jThis++; }   // ONE predictable branch; nothing when off
            if (P<1e-300) P=0.0;          // underflow guard, not a screen
        }
        if (census && P==0.0) { pc.dead++; pc.deadWork+=jThis; }
        sum+=P;
        if (long(i)==target) Pt=P;
    }
    assert(sum>0);
    return Pt/sum;
}

//! \brief THE DERIVATIVE of \c PeriodicBeckeWeight beside it (GeometryOptimizationPlan §13.2; lane B 2026-11-19): ∂w/∂X for every gathered
//! image X in \a im (P-set marked, as the weight left them), the point riding its own centre im[\a target] at the fixed \a offset (r =
//! X_target + offset).  The competitor set and the P-set are HELD (their eps screens are step functions of the geometry: an eps-level term,
//! dropped as every screen's is).  Enzyme over the weight's own product-and-sum (the same arithmetic \c PeriodicBeckeWeight returns).  Sums to
//! zero over the images (translation).  CLIENT: an atom-centred quadrature's weight response (the XC term's grid response on a periodic mesh).
std::vector<rvec3_t> PeriodicBeckeWeightGradient(const std::vector<BeckeImage>& im, size_t target, const rvec3_t& offset, int k);

} // namespace qchem::qcMesh

namespace qchem::qcMesh
{

double BeckePartitionWeight(const std::vector<rvec3_t>& R, size_t owner, const rvec3_t& r, int k)
{
    const size_t natom=R.size();
    if (natom==1) return 1.0;
    std::vector<double> dist(natom), P(natom);
    for (size_t i=0; i<natom; i++) dist[i]=norm(r-R[i]);
    for (size_t i=0; i<natom; i++)
    {
        double Pi=1.0;
        for (size_t j=0; j<natom; j++) if (j!=i)
        {
            double Rij=norm(R[i]-R[j]);
            double mu = (Rij>0.0) ? (dist[i]-dist[j])/Rij : 0.0;  // coincident atoms -> mu=0
            Pi *= BeckeCutoff(mu,k);
        }
        P[i]=Pi;
    }
    double sum=0.0;
    for (size_t i=0; i<natom; i++) sum+=P[i];
    return (sum>0.0) ? P[owner]/sum : 0.0;
}

// THE LEAF'S DERIVATIVE, registered right below it.  Enzyme's entry takes the weight as a function of the positions alone (flattened,
// active), the owner, its offset and the order as a constant context; it differentiates BeckePartitionWeight, BeckeCutoff and norm in
// this unit's IR (doc/Records/EnzymeSpike.md).  Struct and kernel have module linkage, not internal: Enzyme's entry is a template
// declared and never defined, and an internal-linkage argument type makes that ill-formed.
struct BeckeWeightCtx { size_t n, owner; rvec3_t offset; int k; };
double BeckeWeightOfPositions(const double* Rflat, const BeckeWeightCtx* c)
{
    std::vector<rvec3_t> R(c->n);
    for (size_t i=0;i<c->n;i++) R[i]=rvec3_t(Rflat[3*i], Rflat[3*i+1], Rflat[3*i+2]);
    return BeckePartitionWeight(R, c->owner, R[c->owner]+c->offset, c->k);
}

std::vector<rvec3_t> BeckePartitionGradient(const std::vector<rvec3_t>& R, size_t owner, const rvec3_t& r, int k)
{
    const size_t n=R.size();
    std::vector<rvec3_t> out(n, rvec3_t(0,0,0));
    // One atom (w = 1), or the point AT its owner's nucleus: there μ_Aj = -1, where every cutoff is flat (s' = 0), so the gradient is
    // exactly zero -- and |r - R_A| = 0 would put d√x/dx = ∞ into the reverse pass (0·∞ = NaN).
    const rvec3_t offset=r-R[owner];
    if (n==1 || norm(offset)==0.0) return out;
    std::vector<double> x(3*n), g(3*n, 0.0);
    for (size_t i=0;i<n;i++) { x[3*i]=R[i].x; x[3*i+1]=R[i].y; x[3*i+2]=R[i].z; }
    const BeckeWeightCtx c{n, owner, offset, k};
    __enzyme_autodiff<void>(reinterpret_cast<void*>(BeckeWeightOfPositions), enzyme_dup, x.data(), g.data(), enzyme_const, &c);
    for (size_t i=0;i<n;i++) out[i]=rvec3_t(g[3*i], g[3*i+1], g[3*i+2]);
    return out;
}

// THE SCREENED LEAF: the P-set's weight as a function of the image positions alone (flattened, active); the P-set mask, the target and the
// point's offset are constant context.  The same loop as PeriodicBeckeWeight's final double loop (the early exit and underflow guard
// included), so Enzyme differentiates the arithmetic the weight returned.
struct ScreenedBeckeCtx { size_t n, target; const char* inP; rvec3_t offset; int k; };
double ScreenedBeckeWeightOfPositions(const double* Xflat, const ScreenedBeckeCtx* c)
{
    const size_t n=c->n;
    std::vector<rvec3_t> X(n);
    for (size_t i=0;i<n;i++) X[i]=rvec3_t(Xflat[3*i], Xflat[3*i+1], Xflat[3*i+2]);
    const rvec3_t r=X[c->target]+c->offset;
    std::vector<double> d(n);
    for (size_t i=0;i<n;i++) d[i]=norm(r-X[i]);
    double sum=0, Pt=0;
    for (size_t i=0; i<n; i++)
    {
        if (!c->inP[i]) continue;
        double P=1.0;
        for (size_t j=0; j<n && P>0; j++)
        {
            if (j==i) continue;
            double Rij=norm(X[i]-X[j]);
            double mu = (Rij>0.0) ? (d[i]-d[j])/Rij : 0.0;
            P*=BeckeCutoff(mu,c->k);
            if (P<1e-300) P=0.0;
        }
        sum+=P;
        if (i==c->target) Pt=P;
    }
    return (sum>0.0) ? Pt/sum : 0.0;
}

std::vector<rvec3_t> PeriodicBeckeWeightGradient(const std::vector<BeckeImage>& im, size_t target, const rvec3_t& offset, int k)
{
    const size_t n=im.size();
    std::vector<rvec3_t> out(n, rvec3_t(0,0,0));
    // At the owner's nucleus (offset 0) every cutoff is flat in μ = -1 and |r - R| = 0 would put ∞ into the reverse pass: exactly zero.
    if (n<=1 || norm(offset)==0.0) return out;
    if (target>=n) throw std::invalid_argument("PeriodicBeckeWeightGradient: the target is not among the images");
    // THE ACTIVE SET (cost, not a screen of the value: measured 2026-11-19, the AD over every gathered image was half a GPW Becke force).
    // An image j outside the P-set enters only through the factors s(μ_ij), i in the P-set; where every one of them sits on the cutoff's
    // plateau (s = 1 to the last bit, μ_ij <= -1 + δ, the iterated polynomial's flat end) its derivative is zero to that precision, so it is
    // left out of the AD and its factors held at their value, 1.  The P-set and every image that bends some factor stay active.
    const rvec3_t r=im[target].R+offset;
    std::vector<size_t> act;
    act.reserve(n);
    for (size_t j=0;j<n;j++)
    {
        bool keep=im[j].inPSet;
        if (!keep)
        {
            const double dj=norm(r-im[j].R);
            for (size_t i=0;i<n && !keep;i++)
            {
                if (!im[i].inPSet) continue;
                const double Rij=norm(im[i].R-im[j].R);
                const double mu = (Rij>0.0) ? (norm(r-im[i].R)-dj)/Rij : 0.0;
                if (BeckeCutoff(mu,k)!=1.0) keep=true;
            }
        }
        if (keep) act.push_back(j);
    }
    size_t t=act.size();
    for (size_t a=0;a<act.size();a++) if (act[a]==target) t=a;
    if (t==act.size()) throw std::logic_error("PeriodicBeckeWeightGradient: the target is not in the P-set");
    const size_t m=act.size();
    std::vector<double> x(3*m), g(3*m, 0.0);
    std::vector<char> inP(m);
    for (size_t a=0;a<m;a++) { const BeckeImage& b=im[act[a]]; x[3*a]=b.R.x; x[3*a+1]=b.R.y; x[3*a+2]=b.R.z; inP[a]=b.inPSet ? 1 : 0; }
    const ScreenedBeckeCtx c{m, t, inP.data(), offset, k};
    __enzyme_autodiff<void>(reinterpret_cast<void*>(ScreenedBeckeWeightOfPositions), enzyme_dup, x.data(), g.data(), enzyme_const, &c);
    for (size_t a=0;a<m;a++) out[act[a]]=rvec3_t(g[3*a], g[3*a+1], g[3*a+2]);
    return out;
}

} // namespace qchem::qcMesh
