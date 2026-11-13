// File: GaussianRF_AD.C  THE CARTESIAN GAUSSIAN DERIVATIVES, Enzyme's (GeometryOptimizationPlan, "DERIVATIVES ARE ENZYME'S"): the basis
// gradient in r, and the 2-centre S / <p^2> kernels' derivative in the ket centre (KetGradient, below).
//
//
// Enzyme is required (user ruling 2026-11-09), so the hand ∇_r of PrimGaussian, GaussianRF and Polarization is gone: ∇_r[p(r-C) g(r)]
// is the reverse-mode derivative of the value the evaluators already compute, the inline leaves GaussianRF::operator(),
// PrimGaussian::operator() and Polarization::operator() of the module interfaces.  This unit carries the plugin
// (qchem_enzyme_sources in ../../CMakeLists.txt), so the integral kernels in GaussianRF.C keep the tree's flags.
//
// The derivative keeps its finite-difference gates, as every Enzyme derivative does (EnzymeAD/Enzyme#3747, a silently wrong
// reverse-mode gradient, was fixed only at the cdd0d1c2 pin): tests/M_GeometryDerivatives.C, every shell type, near and far, on and off
// the axes through the centre.  #3747's trigger, a callee
// filling a loop-local array through an inner loop with a data-dependent trip count, is absent here: both loops accumulate a scalar.
module;
#include <enzyme/enzyme>
#include <vector>
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


// THE 2-CENTRE KERNELS' KET-CENTRE DERIVATIVE (lane B 2026-11-13; it replaced the libcint seam LibcintPullbackWorkaround, the user's
// 2026-11-10 ruling: "libCint and our own MnD code should not be entangled").  The value is the cache-free twin of GaussianRF::Overlap2C /
// Grad2: per primitive pair the Gaussian product and the Hermite block are BUILT here, from the exponents and the two centres (the cached
// path reads them from the global Ω cache, keyed by primitive ids, which AD would see as constants), and the shared formulas PairOverlap /
// PairGrad2 of the module interface contract them.  So the energy series and its derivative are the same formula, MnD's, and libcint is
// out of the force path.  The Hermite block is the trigger shape of EnzymeAD/Enzyme#3747 (a loop-local array filled through inner loops
// with data-dependent trip counts), fixed at the cdd0d1c2 pin; the derivative keeps its finite-difference gates (tests/M_LatticePullback.C:
// every atom along every axis, all nine virial components, on uncontracted, contracted and d shells).
struct KetPairCtx
{
    GaussianRF::TwoCentre K;
    const double* ea; const double* ca; size_t nprimA; int LA; const double* A;   // the bra radial: exponents, coefficients, L, centre
    const double* eb; const double* cb; size_t nprimB; int LB;                    // the ket radial (its centre is the active argument)
    const Polarization* pa; size_t na; const Polarization* pb; size_t nb;         // the components
    const double* W;                                                              // na x nb weights, row-major
};
double KetPairAt(const double* x, const KetPairCtx* c)
{
    const rvec3_t A(c->A[0],c->A[1],c->A[2]), B(x[0],x[1],x[2]);
    double s=0.0;
    for (size_t i=0;i<c->nprimA;i++)
        for (size_t j=0;j<c->nprimB;j++)
        {
            const GaussProduct gp(c->ea[i],A,c->eb[j],B);
            const Hermite2 H2(gp.p,gp.PA,gp.PB,c->LA+1,c->LB+1);   // as Ω's: one order up for <p^2>
            double t=0.0;
            for (size_t p=0;p<c->na;p++)
                for (size_t q=0;q<c->nb;q++)
                {
                    const double w=c->W[p*c->nb+q];
                    if (w==0.0) continue;
                    t += w * (c->K==GaussianRF::TwoCentre::Overlap ? PairOverlap(gp.p,gp.Eij,H2,c->pa[p],c->pb[q])
                                                                   : PairGrad2(gp.p,gp.Eij,gp.a,gp.b,gp.ab,H2,c->pa[p],c->pb[q]));
                }
            s += c->ca[i]*c->cb[j]*t;
        }
    return s;
}

rvec3_t GaussianRF::KetGradient(TwoCentre K, rf_t& rb, const rvec3_t& B, const Polarization* pa, size_t na, const Polarization* pb,
                                size_t nb, const double* W) const
{
    std::vector<double> ea(itsPrims.size()), eb(rb.itsPrims.size());
    for (size_t i=0;i<ea.size();i++) ea[i]=itsPrims[i]->GetExponent();
    for (size_t j=0;j<eb.size();j++) eb[j]=rb.itsPrims[j]->GetExponent();
    const double A[3]={itsCenter.x,itsCenter.y,itsCenter.z};
    const KetPairCtx c{K, ea.data(), itsCoeff.data(), ea.size(), itsL, A, eb.data(), rb.itsCoeff.data(), eb.size(), rb.itsL,
                       pa, na, pb, nb, W};
    double x[3]={B.x,B.y,B.z}, g[3]={0.0,0.0,0.0};
    __enzyme_autodiff<void>(reinterpret_cast<void*>(KetPairAt), enzyme_dup, x, g, enzyme_const, &c);
    return rvec3_t(g[0],g[1],g[2]);
}

} //namespace qchem::BasisSet::Gaussian::Evaluators::PG_Cart_MnD
