// File: Calculation/GeometryOptimizer.C
//
// qchem::GeometryOptimizer -- a structure-neutral quasi-Newton optimizer that drives a PotentialEnergySurface through a CoordinateMap
// (doc/GeometryOptimizationPlan.md §1; user rulings 2026-11-07, §12: Cartesian behind a coordinate map, the force threshold derived from
// eps by a named calibration, no finite-difference forces).  It sees only q, E(q) and dE/dq.  A molecule maps Cartesian (the external
// modes projected out); a solid will map its cell and free Wyckoff parameters through the same face.
//
// THE METHOD.  BFGS on the Hessian in q, a rational-function (RFO) step in the map's step space, capped by a trust radius that grows on
// a well-predicted step and shrinks on a poor one; a step that RAISES the energy is rejected and the radius cut.  No line search: every
// evaluation is a full SCF point, and the trust radius is the cheaper safeguard.  Tested on analytic model surfaces only until the
// analytic forces exist (ruling §12.3).
module;
#include <cmath>
#include <optional>
#include <string>
#include <vector>
#include <stdexcept>
export module qchem.GeometryOptimizer;
export import qchem.PotentialEnergySurface;   // the surface it drives, PESGeometry, PESFailure
export import qchem.Types;                    // rvec_t, rmat_t, rvec3_t, rmat3d_t
import qchem.Blaze;                           // eigen
import qchem.Reporting;                       // one row per step

export namespace qchem
{

//! \brief How a geometry is parameterised for the optimizer: q <-> (cell, R), forces -> dE/dq, and the subspace a step may move in.
//! Data-free and abstract (CLAUDE.md design rule): the optimizer never sees R, a cell, or what the map does.
class CoordinateMap
{
public:
    virtual ~CoordinateMap() = default;
    virtual size_t Dimension() const = 0;
    virtual rvec_t Start() const = 0;                         //!< q of the reference geometry
    virtual PESGeometry At(const rvec_t& q) const = 0;        //!< the geometry q stands for
    //! dE/dq at \a q from the point's forces (−dE/dR, Ha/bohr, the reference's atom order): −Jᵀ F, J = ∂R/∂q.
    virtual rvec_t Gradient(const rvec_t& q, const std::vector<rvec3_t>& forces) const = 0;
    //! \a v with the directions the map holds fixed at \a q removed: an ORTHOGONAL projection (P = Pᵀ = P²), applied to every gradient,
    //! step and Hessian-update vector, so the optimizer cannot move what the map fixes.
    virtual rvec_t Project(const rvec_t& q, const rvec_t& v) const = 0;
};

//! One point-group operation for \c CartesianMap's symmetrization: atom a goes to atom \c perm[a], at \c rotation · R_a (about the origin).
struct SymmetryOperation
{
    rmat3d_t            rotation;
    std::vector<size_t> perm;
};

//! \brief THE MOLECULAR MAP: q = the Cartesian positions (bohr), flattened atom by atom.
//!
//! \c External::Rigid (a molecule's default) removes the 3 translations and 3 rotations (2 for a linear molecule) about the current
//! centroid, so a step never moves the molecule as a whole; \c None keeps all 3N directions (a model surface that is not
//! translation invariant).  Given \a ops, the step space is also SYMMETRIZED onto the group's totally symmetric subspace,
//! (1/|G|) Σ_g g·v, so a symmetric start stays on the group whatever asymmetric noise its forces carry (level 1, imposed relaxation).
//! The rigid subspace is invariant under the group, so the two projections commute.
class CartesianMap : public virtual CoordinateMap
{
public:
    enum class External {Rigid, None};
    //! THROWS when the reference has no atoms, carries a cell (a solid needs its own map), or is not symmetric under \a ops to \a symTol.
    CartesianMap(const PESGeometry& reference, External external=External::Rigid, std::vector<SymmetryOperation> ops={}, double symTol=1e-8);
    size_t      Dimension() const override {return 3*itsRef.R.size();}
    rvec_t      Start() const override;
    PESGeometry At(const rvec_t& q) const override;
    rvec_t      Gradient(const rvec_t& q, const std::vector<rvec3_t>& forces) const override;
    rvec_t      Project(const rvec_t& q, const rvec_t& v) const override;
    //! The largest |g·R_a − R_{g(a)}| over the operations at \a q (bohr): 0 on the group.  A test and report quantity.
    double      Asymmetry(const rvec_t& q) const;
private:
    PESGeometry                    itsRef;
    External                       itsExternal;
    std::vector<SymmetryOperation> itsOps;
};

//! \brief The named calibration from eps to the force threshold (ruling §12.2).  RIGOROUS: none -- the XC mesh has no error bound.
//! EMPIRICAL: the grid's force error is ~ eta x max|ΔH_ij| (H₂O, LDA/GTH/DZVP-MOLOPT, five rungs, two coordinates: 7-14x; doc/
//! GeometryOptimizationPlan.md §8, ~/Code/qchem6-runs/H2O/fdforces).  ONE molecule: re-measure on analytic forces before trusting it.
struct ForceThresholdCalibration
{
    double eta    = 15.0;   //!< grid force error / max|ΔH_ij| (H₂O: 10, 14, 14, 7 on four rungs)
    double margin = 3.0;    //!< the threshold sits this far above the expected grid force error, so the optimizer never chases the grid
    double ForceThreshold(double eps) const {return margin*eta*eps;}   //!< Ha/bohr: 4.5e-4 at eps = 1e-5
};

//! The stopping thresholds and the trust region.  Converged when max|g| <= maxGradient, rms(g) <= rmsGradient, and either
//! max|step| <= maxStep or |ΔE| <= energyChange.  Units Ha, bohr (q's units for a non-Cartesian map).
struct OptimizerOptions
{
    double maxGradient  = 4.5e-4;
    double rmsGradient  = 3.0e-4;
    double maxStep      = 1.8e-3;
    double energyChange = 1.0e-6;
    double trustInitial = 0.3;
    double trustMax     = 1.0;
    double trustMin     = 1.0e-6;     //!< below it the optimizer stops: \c OptimizeFailure::Why::TrustCollapsed
    double initialCurvature = 0.5;    //!< H₀ = this x 1, Ha/bohr² (until the Lindh model Hessian lands)
    size_t maxSteps     = 200;
    //! The thresholds DERIVED from eps (ruling §12.2): max|g| from the calibration, rms 2/3 of it, the step 4x it (the conventional ratios,
    //! Gaussian's 4.5e-4 / 3.0e-4 / 1.8e-3), and the energy change eps.  The trust region keeps its defaults.
    static OptimizerOptions FromEps(double eps, const ForceThresholdCalibration& cal={});
};

//! One evaluation of the optimizer's history (a VALUE: the report and the tests read it).
struct OptimizationStep
{
    size_t index        = 0;
    double energy       = 0.0;
    double maxGradient  = 0.0;   //!< of the PROJECTED gradient
    double rmsGradient  = 0.0;
    double stepLength   = 0.0;   //!< |s| of the step that led here (0 at the start)
    double trust        = 0.0;   //!< the radius the step was taken under
    bool   accepted     = true;  //!< false: the energy rose, the step was undone and the radius cut
    bool   trustLimited = false; //!< the RFO step was longer than the radius and was scaled to it
    size_t scfIterations= 0;     //!< the point's SCF iterations (PESPoint::iterations)
};

struct OptimizedGeometry
{
    rvec_t                        q;
    PESGeometry                   geometry;
    double                        energy = 0.0;
    std::vector<OptimizationStep> history;   //!< every evaluation in order, the start first
    size_t Steps() const;                    //!< accepted steps
};

struct OptimizeFailure
{
    enum class Why {PES, MaxSteps, TrustCollapsed};
    Why                           why = Why::PES;
    std::string                   details;
    std::optional<PESFailure>     pes;       //!< the point's own failure, when \c why is \c PES
    std::vector<OptimizationStep> history;
};

//! \brief Minimise the surface over the map's coordinates.  Every point is asked for forces; a point that cannot give them (\c NoGradient)
//! or fails its SCF ends the run with that failure named.
class GeometryOptimizer
{
public:
    GeometryOptimizer(PotentialEnergySurface& surface, const CoordinateMap& map, const OptimizerOptions& options={});
    Outcome<OptimizedGeometry, OptimizeFailure> Run();
private:
    PotentialEnergySurface& itsSurface;
    const CoordinateMap&    itsMap;
    OptimizerOptions        itsOpts;
};

} // namespace qchem

namespace qchem
{

namespace
{
double  Comp(const rvec3_t& v, size_t k) {return k==0 ? v.x : k==1 ? v.y : v.z;}
rvec3_t Unit(size_t k)                   {return rvec3_t(k==0, k==1, k==2);}
double  Norm(const rvec_t& v)            {return blazem::norm(v);}
double  Dot(const rvec_t& a, const rvec_t& b) {return blazem::dot(a, b);}
}

//------------------------------------------------------------------------------------------------------------------------------------
// CartesianMap

CartesianMap::CartesianMap(const PESGeometry& reference, External external, std::vector<SymmetryOperation> ops, double symTol)
    : itsRef(reference), itsExternal(external), itsOps(std::move(ops))
{
    if (itsRef.R.empty()) throw std::invalid_argument("CartesianMap: the reference has no atoms");
    if (itsRef.cell)      throw std::invalid_argument("CartesianMap: the reference carries a cell -- a solid needs its own map (cell + Wyckoff)");
    for (size_t k=0;k<itsOps.size();++k)
        if (itsOps[k].perm.size()!=itsRef.R.size())
            throw std::invalid_argument("CartesianMap: operation "+std::to_string(k)+" permutes "+std::to_string(itsOps[k].perm.size())+" atoms, the reference has "
                                        +std::to_string(itsRef.R.size()));
    if (const double a=Asymmetry(Start()); a>symTol)
        throw std::invalid_argument("CartesianMap: the reference is not symmetric under the stated operations (off by "+std::to_string(a)+" bohr)");
}

rvec_t CartesianMap::Start() const
{
    rvec_t q(Dimension());
    for (size_t a=0;a<itsRef.R.size();++a) for (size_t k=0;k<3;++k) q[3*a+k]=Comp(itsRef.R[a], k);
    return q;
}

PESGeometry CartesianMap::At(const rvec_t& q) const
{
    if (q.size()!=Dimension()) throw std::invalid_argument("CartesianMap::At: q has "+std::to_string(q.size())+" entries for "+std::to_string(Dimension()));
    PESGeometry g=itsRef;
    for (size_t a=0;a<g.R.size();++a) g.R[a]=rvec3_t(q[3*a], q[3*a+1], q[3*a+2]);
    return g;
}

rvec_t CartesianMap::Gradient(const rvec_t&, const std::vector<rvec3_t>& forces) const
{
    if (forces.size()!=itsRef.R.size()) throw std::invalid_argument("CartesianMap::Gradient: "+std::to_string(forces.size())+" forces for "
                                                                    +std::to_string(itsRef.R.size())+" atoms");
    rvec_t g(Dimension());
    for (size_t a=0;a<forces.size();++a) for (size_t k=0;k<3;++k) g[3*a+k]=-Comp(forces[a], k);
    return g;
}

double CartesianMap::Asymmetry(const rvec_t& q) const
{
    double worst=0.0;
    const PESGeometry g=At(q);
    for (const auto& op : itsOps)
        for (size_t a=0;a<g.R.size();++a) worst=std::max(worst, norm(op.rotation*g.R[a]-g.R[op.perm[a]]));
    return worst;
}

rvec_t CartesianMap::Project(const rvec_t& q, const rvec_t& v) const
{
    const size_t n=itsRef.R.size();
    rvec_t out=v;
    if (!itsOps.empty())   // the totally symmetric part: (1/|G|) Σ_g g·v, with g moving atom a's displacement to atom perm[a]
    {
        rvec_t s(out.size(), 0.0);
        for (const auto& op : itsOps)
            for (size_t a=0;a<n;++a)
            {
                const rvec3_t d=op.rotation*rvec3_t(out[3*a], out[3*a+1], out[3*a+2]);
                for (size_t k=0;k<3;++k) s[3*op.perm[a]+k]+=Comp(d, k);
            }
        out=s/double(itsOps.size());
    }
    if (itsExternal==External::Rigid && n>1)
    {
        // The rigid modes at q: three translations, three rotations about the centroid; Gram-Schmidt drops the null one of a linear molecule.
        const PESGeometry g=At(q);
        rvec3_t c(0,0,0); for (const auto& r : g.R) c+=r; c/=double(n);
        std::vector<rvec_t> modes;
        for (size_t k=0;k<3;++k)
        {
            rvec_t t(3*n, 0.0); for (size_t a=0;a<n;++a) t[3*a+k]=1.0;
            modes.push_back(t);
            const rvec3_t axis=Unit(k);
            rvec_t r(3*n, 0.0);
            for (size_t a=0;a<n;++a) { const rvec3_t d=Cross(axis, g.R[a]-c); for (size_t j=0;j<3;++j) r[3*a+j]=Comp(d, j); }
            modes.push_back(r);
        }
        std::vector<rvec_t> basis;
        for (rvec_t m : modes)
        {
            for (const auto& b : basis) m-=Dot(b, m)*b;
            const double l=Norm(m);
            if (l>1e-8*std::sqrt(double(n))) basis.push_back(m/l);
        }
        for (const auto& b : basis) out-=Dot(b, out)*b;
    }
    return out;
}

//------------------------------------------------------------------------------------------------------------------------------------
// OptimizerOptions, OptimizedGeometry

OptimizerOptions OptimizerOptions::FromEps(double eps, const ForceThresholdCalibration& cal)
{
    if (!(eps>0.0)) throw std::invalid_argument("OptimizerOptions::FromEps: eps must be > 0");
    OptimizerOptions o;
    o.maxGradient =cal.ForceThreshold(eps);
    o.rmsGradient =o.maxGradient*(2.0/3.0);
    o.maxStep     =o.maxGradient*4.0;
    o.energyChange=eps;
    return o;
}

size_t OptimizedGeometry::Steps() const
{
    size_t n=0;
    for (size_t k=1;k<history.size();++k) if (history[k].accepted) ++n;
    return n;
}

//------------------------------------------------------------------------------------------------------------------------------------
// GeometryOptimizer

GeometryOptimizer::GeometryOptimizer(PotentialEnergySurface& surface, const CoordinateMap& map, const OptimizerOptions& options)
    : itsSurface(surface), itsMap(map), itsOpts(options)
{
    if (!(itsOpts.trustInitial>0.0 && itsOpts.trustMin>0.0 && itsOpts.trustMax>=itsOpts.trustInitial))
        throw std::invalid_argument("GeometryOptimizer: the trust radii must satisfy 0 < trustMin, 0 < trustInitial <= trustMax");
    if (!(itsOpts.initialCurvature>0.0)) throw std::invalid_argument("GeometryOptimizer: the initial curvature must be > 0");
}

namespace
{
double MaxAbs(const rvec_t& v) {double m=0.0; for (double x : v) m=std::max(m, std::fabs(x)); return m;}
double Rms(const rvec_t& v)    {return v.size() ? Norm(v)/std::sqrt(double(v.size())) : 0.0;}

//! An orthonormal basis (the columns) of the map's step space at q: the eigenvectors of the projector with eigenvalue ~ 1.
rmat_t StepBasis(const CoordinateMap& map, const rvec_t& q)
{
    const size_t n=q.size();
    rhmat_t P(n);
    for (size_t j=0;j<n;++j)
    {
        rvec_t e(n, 0.0); e[j]=1.0;
        const rvec_t c=map.Project(q, e);
        for (size_t i=0;i<=j;++i) P(i,j)=c[i];   // P = Pᵀ (an orthogonal projector): the upper triangle carries it
    }
    rvec_t w; rmat_t U;
    blazem::eigen(P, w, U);
    size_t m=0; for (double x : w) if (x>0.5) ++m;
    rmat_t B(n, m);
    size_t k=0;
    for (size_t j=0;j<n;++j) if (w[j]>0.5) { for (size_t i=0;i<n;++i) B(i,k)=U(i,j); ++k; }
    return B;
}
}

Outcome<OptimizedGeometry, OptimizeFailure> GeometryOptimizer::Run()
{
    using O=Outcome<OptimizedGeometry, OptimizeFailure>;
    namespace rpt=qchem::report;
    std::vector<OptimizationStep> history;
    auto evaluate=[&](const rvec_t& q, OptimizationStep& row, double& E, rvec_t& g)->std::optional<PESFailure>
    {
        const PESGeometry geo=itsMap.At(q);
        auto p=itsSurface.ConvergeAt(geo.cell, geo.R, PESRequest{.forces=true});
        if (!p) return p.Error();
        if (!p->forces) return PESFailure{PESFailure::Why::NoGradient, "the surface returned no forces for a forces request", std::nullopt};
        E=p->energy;
        g=itsMap.Project(q, itsMap.Gradient(q, *p->forces));
        row.energy=E; row.maxGradient=MaxAbs(g); row.rmsGradient=Rms(g); row.scfIterations=p->iterations;
        return std::nullopt;
    };
    auto record=[&](const OptimizationStep& s)
    {
        history.push_back(s);
        rpt::Row r("optimizer");
        rpt::Set("step", (long)s.index); rpt::Set("energy", s.energy); rpt::Set("maxGradient", s.maxGradient); rpt::Set("rmsGradient", s.rmsGradient);
        rpt::Set("stepLength", s.stepLength); rpt::Set("trust", s.trust); rpt::Set("accepted", s.accepted); rpt::Set("trustLimited", s.trustLimited);
        rpt::Set("scfIterations", (long)s.scfIterations);
    };
    auto fail=[&](OptimizeFailure::Why why, std::string details, std::optional<PESFailure> pes={})
    {
        return O::Fail(OptimizeFailure{why, std::move(details), std::move(pes), history});
    };
    auto done=[&](const rvec_t& q, double E) {return O::Ok(OptimizedGeometry{q, itsMap.At(q), E, history});};

    const size_t n=itsMap.Dimension();
    rvec_t q=itsMap.Start();
    if (q.size()!=n) throw std::logic_error("GeometryOptimizer: the map's start has "+std::to_string(q.size())+" coordinates for its dimension "+std::to_string(n));
    double E=0.0; rvec_t g;
    OptimizationStep s0; s0.trust=itsOpts.trustInitial;
    if (auto f=evaluate(q, s0, E, g)) { record(s0); return fail(OptimizeFailure::Why::PES, "the starting point: "+f->details, *f); }
    record(s0);
    if (s0.maxGradient<=itsOpts.maxGradient && s0.rmsGradient<=itsOpts.rmsGradient) return done(q, E);   // already a stationary point

    rmat_t H(n, n, 0.0);
    for (size_t i=0;i<n;++i) H(i,i)=itsOpts.initialCurvature;
    double trust=itsOpts.trustInitial;
    for (size_t index=1; index<=itsOpts.maxSteps; ++index)
    {
        // THE RFO STEP in the step space: the lowest eigenvector of [[Hr, gr], [grᵀ, 0]] gives s = v/v_last, a descent direction for any Hr.
        const rmat_t B=StepBasis(itsMap, q);
        const size_t m=B.columns();
        if (m==0) return done(q, E);   // nothing may move
        const rmat_t Hr=blazem::trans(B)*H*B;
        const rvec_t gr=blazem::trans(B)*g;
        rhmat_t A(m+1);
        for (size_t i=0;i<m;++i) { for (size_t j=i;j<m;++j) A(i,j)=0.5*(Hr(i,j)+Hr(j,i)); A(i,m)=gr[i]; }
        A(m,m)=0.0;
        rvec_t w; rmat_t V;
        blazem::eigen(A, w, V);
        rvec_t sr(m);
        if (std::fabs(V(m,0))>1e-12) for (size_t i=0;i<m;++i) sr[i]=V(i,0)/V(m,0);
        else                         sr=-gr;   // the RFO's degenerate case: fall back on steepest descent, which the radius then caps
        OptimizationStep row; row.index=index; row.trust=trust;
        if (const double l=Norm(sr); l>trust) { sr*=trust/l; row.trustLimited=true; }
        rvec_t step=itsMap.Project(q, B*sr);
        const double predicted=Dot(g, step)+0.5*Dot(step, rvec_t(H*step));
        row.stepLength=Norm(step);

        const rvec_t q1=q+step;
        double E1=0.0; rvec_t g1;
        if (auto f=evaluate(q1, row, E1, g1)) { row.accepted=false; record(row); return fail(OptimizeFailure::Why::PES, "step "+std::to_string(index)+": "+f->details, *f); }

        // BFGS on the step that was TAKEN, accepted or not: a rejected step's (s, y) is curvature measured all the same, and dropping it
        // wastes the evaluation (the 6-d bowl from 1.5 bohr: 3 rejections without this update, 2 with it, 2026-11-07).
        // H += y yᵀ/(yᵀs) − (Hs)(Hs)ᵀ/(sᵀHs), skipped when the curvature condition fails, which keeps H positive definite.
        {
            const rvec_t y=itsMap.Project(q1, g1-g);
            const double sy=Dot(step, y);
            if (sy>1e-12*Norm(step)*Norm(y))
            {
                const rvec_t Hs=H*step;
                const double sHs=Dot(step, Hs);
                H+=(y*blazem::trans(y))/sy;
                if (sHs>0.0) H-=(Hs*blazem::trans(Hs))/sHs;
            }
        }
        if (E1>E)   // the energy ROSE: undo the step, cut the radius
        {
            row.accepted=false; record(row);
            trust*=0.25;
            if (trust<itsOpts.trustMin) return fail(OptimizeFailure::Why::TrustCollapsed, "the trust radius fell below "+std::to_string(itsOpts.trustMin)
                                                    +" after an energy rise at step "+std::to_string(index));
            continue;
        }
        record(row);
        const double ratio = predicted<0.0 ? (E1-E)/predicted : 0.0;
        if      (ratio>0.75 && row.trustLimited) trust=std::min(2.0*trust, itsOpts.trustMax);
        else if (ratio<0.25)                     trust=std::max(0.5*trust, itsOpts.trustMin);
        const double dE=E1-E;
        q=q1; E=E1; g=g1;
        const bool forces = row.maxGradient<=itsOpts.maxGradient && row.rmsGradient<=itsOpts.rmsGradient;
        const bool moved  = MaxAbs(step)<=itsOpts.maxStep || std::fabs(dE)<=itsOpts.energyChange;
        if (forces && moved) return done(q, E);
    }
    return fail(OptimizeFailure::Why::MaxSteps, "not converged after "+std::to_string(itsOpts.maxSteps)+" steps");
}

} // namespace qchem
