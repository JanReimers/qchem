// File: GaussianRF_AD.C  THE CARTESIAN GAUSSIAN BASIS GRADIENT, Enzyme's (GeometryOptimizationPlan, "DERIVATIVES ARE ENZYME'S").
//
// Enzyme is required (user ruling 2026-11-09), so the hand ∇_r of PrimGaussian, GaussianRF and Polarization is gone: ∇_r[p(r-C) g(r)]
// is the reverse-mode derivative of the value the evaluators already compute, the inline leaves GaussianRF::operator(),
// PrimGaussian::operator() and Polarization::operator() of the module interfaces.  This unit carries the plugin
// (qchem_enzyme_sources in ../../CMakeLists.txt), so the integral kernels in GaussianRF.C keep the tree's flags.
//
// While EnzymeAD/Enzyme#3747 (a silently wrong reverse-mode gradient) is open the derivative keeps its finite-difference gates
// (tests/M_GeometryDerivatives.C: every shell type, near and far, on and off the axes through the centre).  Its trigger, a callee
// filling a loop-local array through an inner loop with a data-dependent trip count, is absent here: both loops accumulate a scalar.
module;
#include <enzyme/enzyme>
module qchem.BasisSet.Gaussian.Evaluators.PG_Cart_MnD.GaussianRF;

namespace qchem::BasisSet::Gaussian::Evaluators::PG_Cart_MnD
{

// Enzyme's entry takes the value as a function of the point alone (active) and the shell as a constant context.  Struct and kernel have
// module linkage, not internal: the entry is a template declared and never defined, and an internal-linkage argument type makes that
// ill-formed (as in qcMesh's BeckePartition.C).
struct PolarizedRadialCtx { const GaussianRF* rf; Polarization p; };
double PolarizedRadialAt(const double* x, const PolarizedRadialCtx* c)
{
    const rvec3_t r(x[0],x[1],x[2]);
    return c->p(r-c->rf->GetCenter()) * c->rf->GaussianRF::operator()(r);   // qualified: the inline leaf, not the virtual call
}

rvec3_t GaussianRF::Gradient(const Polarization& p, const rvec3_t& r) const
{
    double x[3]={r.x,r.y,r.z}, g[3]={0.0,0.0,0.0};
    const PolarizedRadialCtx c{this, p};
    __enzyme_autodiff<void>(reinterpret_cast<void*>(PolarizedRadialAt), enzyme_dup, x, g, enzyme_const, &c);
    return rvec3_t(g[0],g[1],g[2]);
}

rvec3_t GaussianRF::Gradient(const rvec3_t& r) const
{
    return Gradient(Polarization(0,0,0), r);
}

} //namespace qchem::BasisSet::Gaussian::Evaluators::PG_Cart_MnD
