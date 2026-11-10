// File: Calculation/GeometryOptimizer.C
//
// qchem::GeometryOptimizer -- the GEOMETRY CLIENT of qcOptimize (V-OPTIMIZE-DIP §8 step 2, doc/Records/VOptimizeDIP_Design.md; user
// rulings 2026-11-07, doc/GeometryOptimizationPlan.md §12: Cartesian behind a coordinate map, the force threshold derived from eps by a
// named calibration, no finite-difference forces).  It holds what is GEOMETRY and nothing that is an algorithm:
//   - \c CoordinateMap: q <-> (cell, R), forces -> dE/dq, the step space.  A molecule maps Cartesian (\c CartesianMap, the external
//     modes projected out); a solid maps its free Wyckoff parameters (\c WyckoffMap, qchem.WyckoffMap);
//   - \c PESCost: THE ADAPTER that makes a PotentialEnergySurface seen through a map look like a cost of q (user ruling 2026-11-09,
//     decision 1): R = map.At(q), (E, F) = pes.ConvergeAt(R), g = map.Gradient(q, F);
//   - \c ForceThresholdCalibration and \c OptimizerOptions: the client's stopping rule, derived from eps;
//   - \c GeometryOptimizer: PESCost + the rule + a Stepper (BFGS-RFO by default), run by Optimize::Minimize.
// The algorithm (BFGS + RFO + trust radius, and its noise rules) lives in qcOptimize, which knows no chemistry.
module;
#include <cmath>
#include <memory>
#include <optional>
#include <string>
#include <vector>
#include <stdexcept>
export module qchem.GeometryOptimizer;
export import qchem.PotentialEnergySurface;   // the surface it drives, PESGeometry, PESFailure
export import qchem.Optimize;                 // the cost faces, the stepper, Minimize, BFGS-RFO
export import qchem.Types;                    // rvec_t, rmat_t, rvec3_t, rmat3d_t
import qchem.Blaze;                           // eigen
import qchem.Reporting;                       // the point's SCF iterations, into the minimizer's row

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

//! \brief The geometry client's settings: its STOPPING RULE (converged when max|g| <= maxGradient, rms(g) <= rmsGradient, and either
//! max|step| <= maxStep or |ΔE| <= energyChange; units Ha, bohr -- q's units for a non-Cartesian map), the ENERGY NOISE its points
//! carry, and the default stepper's own knobs.
struct OptimizerOptions
{
    double maxGradient  = 4.5e-4;
    double rmsGradient  = 3.0e-4;
    double maxStep      = 1.8e-3;
    double energyChange = 1.0e-6;
    //! Ha: two evaluations of one geometry may differ by this much (\c PESCost's ValueNoise).  0: exact to roundoff (a model surface).
    double energyNoise  = 0.0;
    size_t maxSteps     = 200;
    Optimize::BFGS_RFO::Options stepper;   //!< the default stepper's trust radii and initial curvature
    //! The thresholds DERIVED from eps (ruling §12.2): max|g| from the calibration, rms 2/3 of it, the step 4x it (the conventional ratios,
    //! Gaussian's 4.5e-4 / 3.0e-4 / 1.8e-3), and the energy change eps.  The energy noise is eps too: the run's grids are sized per point
    //! from eps, so two nearby points' energies agree to about eps and no better.  The stepper keeps its defaults.
    static OptimizerOptions FromEps(double eps, const ForceThresholdCalibration& cal={});
    Optimize::GradientStepValue Rule() const;
};

//! \brief THE ADAPTER: a PotentialEnergySurface seen through a CoordinateMap, as a cost of q.  Every evaluation asks for forces; a
//! point that cannot give them (\c NoGradient) or fails its SCF is the cost's failure, carried back typed (\c PESFailure).  Its
//! capabilities come from what it holds: the step space from the map, the value noise from the client.  Writes the point's SCF
//! iterations into the caller's report row.
class PESCost : public virtual Optimize::tCost<PESFailure>, public virtual Optimize::Projected, public virtual Optimize::ValueNoise
{
public:
    //! THROWS when \a energyNoise < 0.
    PESCost(PotentialEnergySurface& surface, const CoordinateMap& map, double energyNoise);
    size_t Dimension() const override {return itsMap.Dimension();}
    Outcome<Optimize::CostValue, PESFailure> EvaluateTyped(const rvec_t& q, const Optimize::Want&) override;
    rvec_t Project(const rvec_t& q, const rvec_t& v) const override {return itsMap.Project(q, v);}
    double Noise(const rvec_t&, double value) const override;   //!< the larger of the client's energy noise and the value's roundoff
private:
    PotentialEnergySurface& itsSurface;
    const CoordinateMap&    itsMap;
    double                  itsNoise;
};

struct OptimizedGeometry : Optimize::Minimum
{
    PESGeometry geometry;
};
using OptimizeFailure = Optimize::MinimizeFailure<PESFailure>;   //!< \c why Cost carries the point's own \c PESFailure

//! \brief Minimise the surface over the map's coordinates: PESCost + the options' rule + a stepper (BFGS-RFO unless one is given).
class GeometryOptimizer
{
public:
    GeometryOptimizer(PotentialEnergySurface& surface, const CoordinateMap& map, const OptimizerOptions& options={});
    //! With a chosen algorithm: any Stepper minimises the same cost (the LSP claim of V-OPTIMIZE-DIP §8 step 3).
    GeometryOptimizer(PotentialEnergySurface& surface, const CoordinateMap& map, const OptimizerOptions& options,
                      std::unique_ptr<Optimize::Stepper> stepper);
    Outcome<OptimizedGeometry, OptimizeFailure> Run();
private:
    PESCost                            itsCost;
    const CoordinateMap&               itsMap;
    OptimizerOptions                   itsOpts;
    std::unique_ptr<Optimize::Stepper> itsStepper;
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
    if (itsRef.cell)      throw std::invalid_argument("CartesianMap: the reference carries a cell -- a solid maps through WyckoffMap");
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
// OptimizerOptions

OptimizerOptions OptimizerOptions::FromEps(double eps, const ForceThresholdCalibration& cal)
{
    if (!(eps>0.0)) throw std::invalid_argument("OptimizerOptions::FromEps: eps must be > 0");
    OptimizerOptions o;
    o.maxGradient =cal.ForceThreshold(eps);
    o.rmsGradient =o.maxGradient*(2.0/3.0);
    o.maxStep     =o.maxGradient*4.0;
    o.energyChange=eps;
    o.energyNoise =eps;
    return o;
}

Optimize::GradientStepValue OptimizerOptions::Rule() const
{
    return Optimize::GradientStepValue(maxGradient, rmsGradient, maxStep, energyChange, maxSteps);
}

//------------------------------------------------------------------------------------------------------------------------------------
// PESCost

PESCost::PESCost(PotentialEnergySurface& surface, const CoordinateMap& map, double energyNoise)
    : itsSurface(surface), itsMap(map), itsNoise(energyNoise)
{
    if (!(itsNoise>=0.0)) throw std::invalid_argument("PESCost: the energy noise must be >= 0");
}

Outcome<Optimize::CostValue, PESFailure> PESCost::EvaluateTyped(const rvec_t& q, const Optimize::Want&)
{
    using O=Outcome<Optimize::CostValue, PESFailure>;
    const PESGeometry geo=itsMap.At(q);
    auto p=itsSurface.ConvergeAt(geo.cell, geo.R, PESRequest{.forces=true});
    if (!p) return O::Fail(p.Error());
    report::Set("scfIterations", (long)p->iterations);
    if (!p->forces) return O::Fail(PESFailure{PESFailure::Why::NoGradient, "the surface returned no forces for a forces request", std::nullopt});
    return O::Ok(Optimize::CostValue{p->energy, itsMap.Gradient(q, *p->forces)});
}

double PESCost::Noise(const rvec_t&, double value) const {return std::max(itsNoise, Optimize::RoundoffNoise(value));}

//------------------------------------------------------------------------------------------------------------------------------------
// GeometryOptimizer

GeometryOptimizer::GeometryOptimizer(PotentialEnergySurface& surface, const CoordinateMap& map, const OptimizerOptions& options)
    : GeometryOptimizer(surface, map, options, std::make_unique<Optimize::BFGS_RFO>(options.stepper)) {}

GeometryOptimizer::GeometryOptimizer(PotentialEnergySurface& surface, const CoordinateMap& map, const OptimizerOptions& options,
                                     std::unique_ptr<Optimize::Stepper> stepper)
    : itsCost(surface, map, options.energyNoise), itsMap(map), itsOpts(options), itsStepper(std::move(stepper))
{
    if (!itsStepper) throw std::invalid_argument("GeometryOptimizer: no stepper");
}

Outcome<OptimizedGeometry, OptimizeFailure> GeometryOptimizer::Run()
{
    using O=Outcome<OptimizedGeometry, OptimizeFailure>;
    const rvec_t start=itsMap.Start();
    if (start.size()!=itsMap.Dimension())
        throw std::logic_error("GeometryOptimizer: the map's start has "+std::to_string(start.size())+" coordinates for its dimension "+std::to_string(itsMap.Dimension()));
    auto r=Optimize::Minimize(itsCost, *itsStepper, itsOpts.Rule(), start);
    if (!r) return O::Fail(r.Error());
    OptimizedGeometry g;
    static_cast<Optimize::Minimum&>(g)=r.TakeValue();
    g.geometry=itsMap.At(g.q);
    return O::Ok(std::move(g));
}

} // namespace qchem
