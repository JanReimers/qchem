// File: Calculation/WyckoffMap.C
//
// qchem::WyckoffMap -- the level-1 SOLID coordinate map (doc/GeometryOptimizationPlan.md §1; OpenWork §2 "Geometry optimization:
// imposed relaxation"): q are the FREE Wyckoff parameters of the atoms under the imposed space group, the cell held fixed.
// qchem::CellWyckoffMap -- the same, plus the CELL's free strains (the cell half, doc/Records/CellHalf_Design.md §2): see its class.
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
    PESRequest  PointRequest() const override {return PESRequest{.forces=true};}
    rvec_t      Gradient(const rvec_t& q, const PESPoint& point) const override;
    //! Removes the floating-origin modes (rigid translations along a polar axis), which are free parameters but zero modes.
    rvec_t      Project(const rvec_t& q, const rvec_t& v) const override;
    rvec_t      ApplyModelHessian(const rvec_t&, const rvec_t& v, double atomCurvature) const override {return atomCurvature*v;}
    //! −Jᵀ F from the forces alone: the gradient \c Gradient returns, for a caller holding forces but no point (\c CellWyckoffMap).
    rvec_t      ForceGradient(const std::vector<rvec3_t>& forces) const;

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

//! The isotropic elastic model of a cell's curvature, for \c CellWyckoffMap's model Hessian.  A GUESS, stated: it only seeds the
//! quasi-Newton history (BFGS learns the real one), so its job is the right ORDER of magnitude, cell block against atom block.
struct ElasticModel
{
    static constexpr double kHartreePerBohr3InGPa = 29421.02648438959;   //!< 1 Ha/bohr³ in GPa (CODATA 2018)
    double bulkModulus = 100.0/kHartreePerBohr3InGPa;   //!< B, Ha/bohr³: 100 GPa, a typical oxide's (rutile TiO₂ ~210, Si ~98, NaF ~46)
    double poisson     = 0.25;                         //!< ν: λ = 3Bν/(1+ν), μ = 3B(1−2ν)/(2(1+ν)); 0.25 gives λ = μ = 0.6 B
};

//! \brief THE SOLID MAP WITH A CELL (level 1, the cell half; doc/Records/CellHalf_Design.md §2, ruled 2026-11-10): q = (the free
//! Wyckoff parameters, the cell's free strains).
//!
//! THE PARAMETERISATION.  F = 1 + ε is the deformation from the REFERENCE cell, A = F A_ref, and the atoms ride it:
//!     R(q) = F (R_ref + J (q_int − q_int₀)),
//! J the \c WyckoffMap's orthonormal Wyckoff modes at the reference metric.  The fractional coordinates depend on q_int alone, so at
//! fixed q_int a cell coordinate's gradient is the stress alone.  ε is SYMMETRIC (a rotation is not a coordinate) and
//! GROUP-INVARIANT, ε = (1/|G|) Σ_g C_g ε C_gᵀ with C_g the ops' Cartesian linear parts: then F commutes with every C_g, so the
//! strained crystal keeps every op (in fractional coordinates the ops do not change).  No Bravais table: the invariant subspace IS the
//! lattice type -- cubic 1, tetragonal and hexagonal 2, orthorhombic 3, monoclinic 4, triclinic 6.  E_k is its Frobenius-orthonormal
//! basis and ε = Σ_k e_k E_k.
//!
//! UNITS (ruling 3): every coordinate is in bohr.  q_cell,k = L e_k with L = Ω_ref^{1/3}, so a cell step reads like an atom step.  The
//! gradients, from ∂E/∂F = Ω σ F^{−T} (fractional coordinates held; σ the point's stress, \c PESPoint::stress):
//!     g_int     = −Jᵀ Fᵀ F_atoms,
//!     g_cell,k  = (Ω/L) tr(F^{−1} σ E_k),
//! both Ha/bohr.  (The design note's tr(σ F^{−T} E_k) differs off the commuting cases; this is the contraction of ∂E/∂F with
//! ∂F/∂q_k = E_k / L, and its FD gate is the triclinic case.)
//!
//! THE MODEL HESSIAN (ruling 3): the atom block is the client's curvature; the cell block is the isotropic elastic energy
//! ½ Ω_ref [λ (tr ε)² + 2 μ ε:ε], in q: (Ω_ref / L²) [λ t tᵀ + 2 μ 1], t_k = tr E_k; no atom-cell coupling.
//!
//! The floating origin is held as \c WyckoffMap holds it (a polar axis is invariant under every C_g, so F keeps it one).
class CellWyckoffMap : public virtual CoordinateMap
{
public:
    //! THROWS as \c WyckoffMap does.
    CellWyckoffMap(const PESGeometry& reference, const std::vector<Symmetry::SymOp>& ops, const ElasticModel& model={}, double symTol=1e-6);
    CellWyckoffMap(const PESGeometry& reference, const Symmetry::Lattice_3D::SpaceGroup& group, const ElasticModel& model={}, double symTol=1e-6);

    size_t      Dimension() const override {return InternalDimension()+CellDimension();}
    rvec_t      Start() const override;
    PESGeometry At(const rvec_t& q) const override;
    PESRequest  PointRequest() const override {return PESRequest{.forces=true, .stress=true};}
    rvec_t      Gradient(const rvec_t& q, const PESPoint& point) const override;
    rvec_t      Project(const rvec_t& q, const rvec_t& v) const override;
    rvec_t      ApplyModelHessian(const rvec_t& q, const rvec_t& v, double atomCurvature) const override;

    size_t   InternalDimension() const {return itsInternal.Dimension();}   //!< the free Wyckoff parameters (q's first block)
    size_t   CellDimension() const {return itsE.size();}                   //!< the free strains (q's last block): 1..6
    size_t   HeldModes() const {return itsInternal.HeldModes();}
    //! The strain ε (dimensionless, symmetric) that \a q stands for.
    rmat3d_t Strain(const rvec_t& q) const;
    //! The free strain directions E_k (Frobenius-orthonormal, symmetric, group-invariant).
    const std::vector<rmat3d_t>& StrainBasis() const {return itsE;}

private:
    WyckoffMap            itsInternal;
    rmat3d_t              itsAref;
    std::vector<rmat3d_t> itsE;
    double                itsL;
    double                itsOmegaRef;
    ElasticModel          itsModel;
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

rvec_t WyckoffMap::Gradient(const rvec_t&, const PESPoint& point) const
{
    if (!point.forces) throw std::invalid_argument("WyckoffMap::Gradient: the point carries no forces");
    return ForceGradient(*point.forces);
}

rvec_t WyckoffMap::ForceGradient(const std::vector<rvec3_t>& forces) const
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

//------------------------------------------------------------------------------------------------------------------------------------
// CellWyckoffMap

namespace
{
rmat3d_t Zero3() {return rmat3d_t(0,0,0, 0,0,0, 0,0,0);}
rmat3d_t One3()  {return rmat3d_t(1,0,0, 0,1,0, 0,0,1);}
//! The Frobenius-orthonormal basis of the symmetric 3x3 matrices: e_i e_iᵀ, and (e_i e_jᵀ + e_j e_iᵀ)/√2.
std::vector<rmat3d_t> SymmetricBasis()
{
    std::vector<rmat3d_t> B;
    for (size_t i=1;i<=3;++i) { rmat3d_t m=Zero3(); m(i,i)=1.0; B.push_back(m); }
    for (size_t i=1;i<=3;++i) for (size_t j=i+1;j<=3;++j) { rmat3d_t m=Zero3(); m(i,j)=m(j,i)=1.0/std::sqrt(2.0); B.push_back(m); }
    return B;
}
double Frobenius(const rmat3d_t& a, const rmat3d_t& b) {double s=0; for (size_t i=1;i<=3;++i) for (size_t j=1;j<=3;++j) s+=a(i,j)*b(i,j); return s;}
double Trace(const rmat3d_t& a) {return a(1,1)+a(2,2)+a(3,3);}
}

CellWyckoffMap::CellWyckoffMap(const PESGeometry& reference, const Symmetry::Lattice_3D::SpaceGroup& group, const ElasticModel& model, double symTol)
    : CellWyckoffMap(reference, OpsOf(group), model, symTol) {}

CellWyckoffMap::CellWyckoffMap(const PESGeometry& reference, const std::vector<Symmetry::SymOp>& ops, const ElasticModel& model, double symTol)
    : itsInternal(reference, ops, symTol)   // THROWS on every refusal the positions-only map makes (no cell among them)
    , itsAref(*reference.cell), itsModel(model)
{
    if (!(itsModel.bulkModulus>0.0) || !(itsModel.poisson>-1.0 && itsModel.poisson<0.5))
        throw std::invalid_argument("CellWyckoffMap: the elastic model needs B > 0 and −1 < ν < ½");
    itsOmegaRef=std::fabs(Determinant(itsAref));
    itsL=std::cbrt(itsOmegaRef);
    // THE INVARIANT STRAINS: the average of ε → C ε Cᵀ over the group, a symmetric projector on the 6-d symmetric matrices (C orthogonal:
    // the positions-only map checked that the ops preserve the metric).  Its eigenvalue-1 eigenvectors are the free strains.
    const rmat3d_t Ainv=Invert(itsAref);
    const std::vector<rmat3d_t> B=SymmetricBasis();
    rhmat_t P(6);
    for (size_t m=0;m<6;++m) for (size_t n=m;n<6;++n) P(m,n)=0.0;
    for (const auto& op : ops)
    {
        const rmat3d_t C=itsAref*op.W*Ainv;
        for (size_t n=0;n<6;++n)
        {
            const rmat3d_t CBC=C*B[n]*Transpose(C);
            for (size_t m=0;m<=n;++m) P(m,n)+=Frobenius(B[m], CBC)/double(ops.size());
        }
    }
    rvec_t w; rmat_t U;
    blazem::eigen(P, w, U);
    for (size_t e=0;e<6;++e)
    {
        if (w[e]<0.5) continue;
        rmat3d_t E=Zero3();
        for (size_t m=0;m<6;++m) E+=U(m,e)*B[m];
        itsE.push_back(E);
    }
    if (itsE.empty()) throw std::logic_error("CellWyckoffMap: no invariant strain -- the hydrostatic one is invariant under every orthogonal op");
}

rmat3d_t CellWyckoffMap::Strain(const rvec_t& q) const
{
    if (q.size()!=Dimension()) throw std::invalid_argument("CellWyckoffMap: q has "+std::to_string(q.size())+" entries for "+std::to_string(Dimension()));
    rmat3d_t eps=Zero3();
    for (size_t k=0;k<itsE.size();++k) eps+=(q[InternalDimension()+k]/itsL)*itsE[k];
    return eps;
}

rvec_t CellWyckoffMap::Start() const
{
    rvec_t q(Dimension(), 0.0);
    const rvec_t qi=itsInternal.Start();
    for (size_t i=0;i<qi.size();++i) q[i]=qi[i];
    return q;
}

namespace
{
rvec_t InternalPart(const rvec_t& q, size_t n) {rvec_t qi(n); for (size_t i=0;i<n;++i) qi[i]=q[i]; return qi;}
}

PESGeometry CellWyckoffMap::At(const rvec_t& q) const
{
    const rmat3d_t F=One3()+Strain(q);   // THROWS on a wrong size
    PESGeometry g=itsInternal.At(InternalPart(q, InternalDimension()));
    for (auto& r : g.R) r=F*r;
    g.cell=F*itsAref;
    return g;
}

rvec_t CellWyckoffMap::Gradient(const rvec_t& q, const PESPoint& point) const
{
    if (!point.forces) throw std::invalid_argument("CellWyckoffMap::Gradient: the point carries no forces");
    if (!point.stress) throw std::invalid_argument("CellWyckoffMap::Gradient: the point carries no stress");
    const rmat3d_t F=One3()+Strain(q), Ft=Transpose(F), Finv=Invert(F);
    std::vector<rvec3_t> f=*point.forces;
    for (auto& x : f) x=Ft*x;                      // −∂E/∂R_ref-metric = Fᵀ F_atoms
    const rvec_t gi=itsInternal.ForceGradient(f);   // THROWS on a wrong count
    rvec_t g(Dimension(), 0.0);
    for (size_t i=0;i<gi.size();++i) g[i]=gi[i];
    const double Omega=std::fabs(Determinant(F*itsAref));
    const rmat3d_t FinvSigma=Finv*(*point.stress);
    for (size_t k=0;k<itsE.size();++k) g[InternalDimension()+k]=(Omega/itsL)*Trace(FinvSigma*itsE[k]);
    return g;
}

rvec_t CellWyckoffMap::Project(const rvec_t& q, const rvec_t& v) const
{
    const size_t n=InternalDimension();
    const rvec_t pi=itsInternal.Project(InternalPart(q, n), InternalPart(v, n));
    rvec_t out=v;
    for (size_t i=0;i<n;++i) out[i]=pi[i];
    return out;
}

rvec_t CellWyckoffMap::ApplyModelHessian(const rvec_t& q, const rvec_t& v, double atomCurvature) const
{
    if (v.size()!=Dimension() || q.size()!=Dimension()) throw std::invalid_argument("CellWyckoffMap::ApplyModelHessian: wrong size");
    const size_t n=InternalDimension();
    rvec_t out(Dimension());
    for (size_t i=0;i<n;++i) out[i]=atomCurvature*v[i];
    const double B=itsModel.bulkModulus, nu=itsModel.poisson;
    const double lambda=3.0*B*nu/(1.0+nu), mu=3.0*B*(1.0-2.0*nu)/(2.0*(1.0+nu));
    const double scale=itsOmegaRef/(itsL*itsL);
    double tv=0.0;
    for (size_t k=0;k<itsE.size();++k) tv+=Trace(itsE[k])*v[n+k];
    for (size_t k=0;k<itsE.size();++k) out[n+k]=scale*(lambda*Trace(itsE[k])*tv+2.0*mu*v[n+k]);
    return out;
}

} // namespace qchem
