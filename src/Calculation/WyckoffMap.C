// File: Calculation/WyckoffMap.C
//
// qchem::WyckoffMap -- the level-1 SOLID coordinate map (doc/GeometryOptimizationPlan.md §1; OpenWork §2 "Geometry optimization:
// imposed relaxation"): q are the FREE Wyckoff parameters of the atoms under the imposed space group, the cell held fixed (the cell
// half needs stress, the variable-cell seam, and comes later).
//
// THE MAP.  Each orbit of the group (a Wyckoff set) has a representative atom a.  Its SITE GROUP S_a (the ops that fix it modulo the
// lattice, as \c SpaceGroup::SiteStabilizer) leaves a subspace of displacements fixed: the range of P_a = (1/|S_a|) Σ_{s∈S_a} C_s,
// with C = A W A⁻¹ the op's Cartesian linear part (orthogonal, so P_a is an orthogonal projector).  Its orthonormal basis u_1..u_k
// are the representative's free directions -- the Wyckoff letter's x, y, z, up to a rotation within the subspace.  A member b = g·a
// of the orbit then moves by C_g δR_a (the translation τ_g and the lattice vector that brought b home are fixed by the reference),
// so mode i of the orbit is the 3N vector
//     J[:, i]_b = C_{g_b} u_i / √m        (m = the orbit's multiplicity)
// and R(q) = R₀ + J (q − q₀).  The columns of J are ORTHONORMAL (C_g is orthogonal; orbits are disjoint), so q is an isometric
// coordinate of the symmetric subspace: steps, trust radii and gradient norms are Cartesian bohr and Ha/bohr, q = Jᵀ R, and the
// Hessian in q is Jᵀ H_R J.  dE/dq = −Jᵀ F is the forces' totally symmetric part: asymmetric force noise drops out, and R(q) is on
// the group for EVERY q to roundoff, by construction.
//
// THE FLOATING ORIGIN.  In a polar group a rigid translation along the polar axis is a free Wyckoff parameter of every orbit, and a
// zero mode of any periodic energy.  \c Dimension() counts it (it IS a free parameter); \c Project holds it fixed, so the optimizer
// steps in Dimension() minus the polar-axis count.
module;
#include <cmath>
#include <string>
#include <vector>
#include <stdexcept>
export module qchem.WyckoffMap;
export import qchem.GeometryOptimizer;                     // CoordinateMap, PESGeometry
export import qchem.Symmetry.Lattice_3D.SpaceGroup;        // SpaceGroup, SymOp (the imposed group's {W|τ})
import qchem.Blaze;                                        // eigen, rmat_t

export namespace qchem
{

//! \brief THE SOLID MAP, level 1: q = the free Wyckoff parameters of the atoms under the imposed group; the cell is the reference's.
//!
//! The group is handed in as its DIRECT-frame operations \f$\{W|\tau\}\f$ (fractional, \f$f \to W f + \tau\f$): a plain space
//! group's \c Ops(), or a magnetic run's \c ShubnikovOps (the spin action \f$\sigma\f$ acts on field channels, never on positions,
//! so it is ignored here).  It must be the group the SESSION imposes: the factory that builds the session hands both the same ops.
class WyckoffMap : public virtual CoordinateMap
{
public:
    //! THROWS (broken invariants, no caller can act on them): the reference has no cell or no atoms; \a ops is empty or not a closed
    //! group; or the reference is not symmetric under \a ops to \a symTol (bohr) -- an op maps an atom onto no atom of its species.
    WyckoffMap(const PESGeometry& reference, const std::vector<Symmetry::SymOp>& ops, double symTol=1e-6);
    //! The same, from a detected \c SpaceGroup (its \c Ops()).
    WyckoffMap(const PESGeometry& reference, const Symmetry::Lattice_3D::SpaceGroup& group, double symTol=1e-6);

    size_t      Dimension() const override {return itsJ.columns();}   //!< the number of free Wyckoff parameters (floating origin included)
    rvec_t      Start() const override {return itsQ0;}                 //!< q₀ = Jᵀ R₀
    PESGeometry At(const rvec_t& q) const override;
    rvec_t      Gradient(const rvec_t& q, const std::vector<rvec3_t>& forces) const override;
    //! Removes the floating-origin modes (rigid translations along a polar axis), which are free parameters but zero modes.
    rvec_t      Project(const rvec_t& q, const rvec_t& v) const override;

    //! The largest |W f_a + τ − f_{g(a)}| (Cartesian, minimum image, bohr) over the ops and atoms at \a q: 0 on the group to roundoff.
    double Asymmetry(const rvec_t& q) const;
    size_t Orbits() const {return itsOrbitCount;}   //!< the number of Wyckoff sets (orbits of the group on the atoms)
    size_t HeldModes() const {return itsHeld.size();} //!< the floating-origin modes \c Project removes (0 unless the group is polar)

private:
    PESGeometry                     itsRef;
    rmat3d_t                        itsA, itsAinv;
    std::vector<Symmetry::SymOp>    itsOps;
    std::vector<std::vector<size_t>> itsPerm;   //!< itsPerm[g][a] = the atom op g maps atom a onto (fixed at the reference)
    rmat_t                          itsJ;       //!< ∂R/∂q, 3N x D, orthonormal columns
    rvec_t                          itsQ0;
    std::vector<rvec_t>             itsHeld;    //!< orthonormal q-space floating-origin modes
    size_t                          itsOrbitCount=0;
};

} // namespace qchem

namespace qchem
{

namespace
{
double& Comp(rvec3_t& v, size_t k)       {return k==0 ? v.x : k==1 ? v.y : v.z;}
double  Comp(const rvec3_t& v, size_t k) {return k==0 ? v.x : k==1 ? v.y : v.z;}
//! The minimum-image fractional difference (each component reduced to [-½, ½]).
rvec3_t Wrap(rvec3_t d) {for (size_t k=0;k<3;++k) Comp(d,k)-=std::round(Comp(d,k)); return d;}

std::vector<Symmetry::SymOp> OpsOf(const Symmetry::Lattice_3D::SpaceGroup& g)
{
    std::vector<Symmetry::SymOp> ops;
    for (const auto& op : g.Ops()) ops.push_back({op.W, op.tau});
    return ops;
}
}

WyckoffMap::WyckoffMap(const PESGeometry& reference, const Symmetry::Lattice_3D::SpaceGroup& group, double symTol)
    : WyckoffMap(reference, OpsOf(group), symTol) {}

WyckoffMap::WyckoffMap(const PESGeometry& reference, const std::vector<Symmetry::SymOp>& ops, double symTol)
    : itsRef(reference), itsOps(ops)
{
    if (!itsRef.cell)     throw std::invalid_argument("WyckoffMap: the reference has no cell -- a molecule maps through CartesianMap");
    if (itsRef.R.empty()) throw std::invalid_argument("WyckoffMap: the reference has no atoms");
    if (itsOps.empty())   throw std::invalid_argument("WyckoffMap: no operations -- hand it the imposed group (the identity at least)");
    if (itsRef.Z.size()!=itsRef.R.size())
        throw std::invalid_argument("WyckoffMap: "+std::to_string(itsRef.Z.size())+" species for "+std::to_string(itsRef.R.size())+" atoms");
    Symmetry::Lattice_3D::RequireClosedGroup(itsOps, "WyckoffMap");
    itsA=*itsRef.cell;
    itsAinv=Invert(itsA);
    const size_t N=itsRef.R.size();

    // THE PERMUTATION each op makes of the reference's atoms: the nearest atom of the same species to W f_a + τ (minimum image).
    std::vector<rvec3_t> f(N);
    for (size_t a=0;a<N;++a) f[a]=itsAinv*itsRef.R[a];
    itsPerm.assign(itsOps.size(), std::vector<size_t>(N));
    for (size_t g=0;g<itsOps.size();++g)
        for (size_t a=0;a<N;++a)
        {
            const rvec3_t img=itsOps[g].W*f[a]+itsOps[g].tau;
            double best=1e300; size_t at=N;
            for (size_t b=0;b<N;++b)
            {
                if (itsRef.Z[b]!=itsRef.Z[a]) continue;
                const double d=norm(itsA*Wrap(img-f[b]));
                if (d<best) {best=d; at=b;}
            }
            if (at==N || best>symTol)
                throw std::invalid_argument("WyckoffMap: the reference is not symmetric under operation "+std::to_string(g)+": atom "+std::to_string(a)
                                            +" (Z="+std::to_string(itsRef.Z[a])+") lands "+(at==N ? std::string("on no atom of its species")
                                            : std::to_string(best)+" bohr from the nearest one")+" (tolerance "+std::to_string(symTol)+" bohr)");
            itsPerm[g][a]=at;
        }

    // THE ORBITS (Wyckoff sets) and, per member, an op that brings the representative there.
    std::vector<rmat3d_t> C(itsOps.size());
    for (size_t g=0;g<itsOps.size();++g) C[g]=itsA*itsOps[g].W*itsAinv;   // Cartesian linear part: orthogonal for a symmetry of the cell
    std::vector<rvec_t> columns;
    std::vector<bool> seen(N, false);
    for (size_t a=0;a<N;++a)
    {
        if (seen[a]) continue;
        ++itsOrbitCount;
        std::vector<std::pair<size_t,size_t>> members;   // (atom, op taking a to it)
        for (size_t g=0;g<itsOps.size();++g)
            if (const size_t b=itsPerm[g][a]; !seen[b]) {seen[b]=true; members.push_back({b, g});}
        // THE SITE GROUP's fixed subspace: P = (1/|S|) Σ_{s∈S} C_s over the ops that fix a (W f + τ ≡ f mod the lattice; the same set
        // as SpaceGroup::SiteStabilizer).  Symmetric for orthogonal C_s; its eigenvalue-1 eigenvectors are the free directions.
        rhmat_t P(3);
        for (size_t i=0;i<3;++i) for (size_t j=i;j<3;++j) P(i,j)=0.0;
        size_t nS=0;
        for (size_t g=0;g<itsOps.size();++g)
            if (itsPerm[g][a]==a)
            {
                ++nS;
                for (size_t i=0;i<3;++i) for (size_t j=i;j<3;++j) P(i,j)+=0.5*(C[g](i+1,j+1)+C[g](j+1,i+1));
            }
        for (size_t i=0;i<3;++i) for (size_t j=i;j<3;++j) P(i,j)/=double(nS);
        rvec_t w; rmat_t U;
        blazem::eigen(P, w, U);
        const double scale=1.0/std::sqrt(double(members.size()));
        for (size_t e=0;e<3;++e)
        {
            if (w[e]<0.5) continue;
            const rvec3_t u(U(0,e), U(1,e), U(2,e));
            rvec_t col(3*N, 0.0);
            for (const auto& [b, g] : members) { const rvec3_t d=C[g]*u; for (size_t k=0;k<3;++k) col[3*b+k]=scale*Comp(d,k); }
            columns.push_back(col);
        }
    }
    itsJ=rmat_t(3*N, columns.size());
    for (size_t c=0;c<columns.size();++c) for (size_t i=0;i<3*N;++i) itsJ(i,c)=columns[c][i];
    {
        // J's columns are orthonormal by construction; a cell the ops do not preserve (C not orthogonal) would break it -- check once.
        const rmat_t JtJ=blazem::trans(itsJ)*itsJ;
        for (size_t i=0;i<JtJ.rows();++i) for (size_t j=0;j<JtJ.columns();++j)
            if (std::fabs(JtJ(i,j)-(i==j ? 1.0 : 0.0))>1e-8)
                throw std::logic_error("WyckoffMap: the free modes are not orthonormal (JᵀJ off by "+std::to_string(JtJ(i,j)-(i==j))
                                       +"): the operations do not preserve the cell's metric");
    }
    rvec_t R0(3*N);
    for (size_t a=0;a<N;++a) for (size_t k=0;k<3;++k) R0[3*a+k]=Comp(itsRef.R[a], k);
    itsQ0=blazem::trans(itsJ)*R0;

    // THE FLOATING ORIGIN: a rigid translation's component inside the step space (nonzero only along a polar axis).
    for (size_t k=0;k<3;++k)
    {
        rvec_t T(3*N, 0.0); for (size_t a=0;a<N;++a) T[3*a+k]=1.0/std::sqrt(double(N));
        rvec_t t=blazem::trans(itsJ)*T;
        for (const auto& h : itsHeld) t-=blazem::dot(h, t)*h;
        if (const double l=blazem::norm(t); l>1e-6) itsHeld.push_back(t/l);
    }
}

PESGeometry WyckoffMap::At(const rvec_t& q) const
{
    if (q.size()!=Dimension()) throw std::invalid_argument("WyckoffMap::At: q has "+std::to_string(q.size())+" entries for "+std::to_string(Dimension()));
    const rvec_t dR=itsJ*(q-itsQ0);
    PESGeometry g=itsRef;
    for (size_t a=0;a<g.R.size();++a) for (size_t k=0;k<3;++k) Comp(g.R[a], k)+=dR[3*a+k];
    return g;
}

rvec_t WyckoffMap::Gradient(const rvec_t&, const std::vector<rvec3_t>& forces) const
{
    if (forces.size()!=itsRef.R.size()) throw std::invalid_argument("WyckoffMap::Gradient: "+std::to_string(forces.size())+" forces for "
                                                                    +std::to_string(itsRef.R.size())+" atoms");
    rvec_t F(3*forces.size());
    for (size_t a=0;a<forces.size();++a) for (size_t k=0;k<3;++k) F[3*a+k]=Comp(forces[a], k);
    return -(blazem::trans(itsJ)*F);
}

rvec_t WyckoffMap::Project(const rvec_t&, const rvec_t& v) const
{
    rvec_t out=v;
    for (const auto& h : itsHeld) out-=blazem::dot(h, out)*h;
    return out;
}

double WyckoffMap::Asymmetry(const rvec_t& q) const
{
    const PESGeometry g=At(q);
    std::vector<rvec3_t> f(g.R.size());
    for (size_t a=0;a<f.size();++a) f[a]=itsAinv*g.R[a];
    double worst=0.0;
    for (size_t o=0;o<itsOps.size();++o)
        for (size_t a=0;a<f.size();++a) worst=std::max(worst, norm(itsA*Wrap(itsOps[o].W*f[a]+itsOps[o].tau-f[itsPerm[o][a]])));
    return worst;
}

} // namespace qchem
