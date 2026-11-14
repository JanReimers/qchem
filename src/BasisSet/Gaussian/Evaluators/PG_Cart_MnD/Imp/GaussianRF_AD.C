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
#include <stdexcept>
#include <string>
#include <vector>
module qchem.BasisSet.Gaussian.Evaluators.PG_Cart_MnD.GaussianRF;
import qchem.BasisSet.Gaussian.Evaluators.Internal.Boys.AuxillaryFJ;   // F_j (the nuclear twin)
import qchem.BasisSet.Gaussian.Evaluators.Internal.MnD.RNLM;          // HermiteCoulomb
import qchem.BasisSet.Gaussian.Evaluators.Internal.MnD.Triangle3D;    // the layer/row index of a Hermite-Coulomb triangle

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


// THE NUCLEAR ATTRACTION'S GEOMETRY DERIVATIVE, every centre moving (lane B 2026-11-13; it replaced libcint's int1e_ipnuc + int1e_iprinv
// for the molecular PG_Cart pullback).  The cache-free twin of PrimGaussian::Nuclear: per primitive pair the Gaussian product and the
// Hermite2 block, per nucleus F_j(T) (Boys::AuxillaryFJ::GetFjAt) and the Hermite-Coulomb recursion (MnD::HermiteCoulomb) on a flat local
// store, summed with -Z_C as the cached RNLM::Add does, contracted by the shared PairNuclear.  The active input is every centre at once,
// x = [A, B, C_0 .. C_{n-1}], so one reverse pass gives the bra, the ket and every nucleus.  Gated by tests/M_GeometryDerivatives.C
// (the FD of the energy's own matrices; libcint's NuclearGradient as an oracle) and tests/M_KetGradient.C (the kernel's own FD).
struct NucPairCtx
{
    const double* ea; const double* ca; size_t nprimA; int LA;
    const double* eb; const double* cb; size_t nprimB; int LB;
    const Polarization* pa; size_t na; const Polarization* pb; size_t nb;
    const double* W; const double* Z; size_t nC;
    const std::vector<Polarization>* NLMs; const Internal::Boys::AuxillaryFJ* Fj;
};
double NucPairAt(const double* x, const NucPairCtx* c)
{
    using Internal::MnD::Triangle3D;
    const rvec3_t A(x[0],x[1],x[2]), B(x[3],x[4],x[5]);
    const int Ltot=c->LA+c->LB;
    auto tri=[](int n){return size_t(n+1)*(n+2)*(n+3)/6;};
    std::vector<size_t> off(Ltot+2,0);
    for (int j=0;j<=Ltot;j++) off[j+1]=off[j]+tri(Ltot-j);
    std::vector<double> store(off[Ltot+1],0.0), Rsum(tri(Ltot),0.0), Fj(Ltot+1,0.0);
    auto at=[&](int j, int N, int L, int M) -> double&
    {   const int n=Ltot-j; return store[off[j]+Triangle3D::Layer(n-N)+Triangle3D::Row(n-N-L)+size_t(n-N-L-M)]; };
    double s=0.0;
    for (size_t i=0;i<c->nprimA;i++)
        for (size_t j=0;j<c->nprimB;j++)
        {
            const GaussProduct gp(c->ea[i],A,c->eb[j],B);
            const Hermite2 H2(gp.p,gp.PA,gp.PB,c->LA+1,c->LB+1);   // as Ω's
            for (double& r : Rsum) r=0.0;
            for (size_t k=0;k<c->nC;k++)
            {
                const rvec3_t dR=gp.P-rvec3_t(x[6+3*k],x[7+3*k],x[8+3*k]);
                c->Fj->GetFjAt(gp.p*dR*dR, Fj);
                Internal::MnD::HermiteCoulomb(Ltot, gp.p, dR, Fj, at);
                for (size_t q=0;q<Rsum.size();q++) Rsum[q]+=store[q]*(-1.0*c->Z[k]);   // layer j = 0 is store[0 .. tri(Ltot))
            }
            auto Rof=[&](const Polarization& m){ return Rsum[Triangle3D::Layer(Ltot-m.n)+Triangle3D::Row(Ltot-m.n-m.l)+size_t(Ltot-m.n-m.l-m.m)]; };
            double t=0.0;
            for (size_t p=0;p<c->na;p++)
                for (size_t q=0;q<c->nb;q++)
                {
                    const double w=c->W[p*c->nb+q];
                    if (w==0.0) continue;   // a constant weight: safe
                    t += w*PairNuclear(gp.p,gp.Eij,H2,c->pa[p],c->pb[q],*c->NLMs,Rof);
                }
            s += c->ca[i]*c->cb[j]*t;
        }
    return s;
}

void GaussianRF::NuclearGradient(rf_t& rb, const Polarization* pa, size_t na, const Polarization* pb, size_t nb, const double* W,
                                 const std::vector<rvec3_t>& C, const std::vector<double>& Z, rvec3_t& gA, rvec3_t& gB,
                                 std::vector<rvec3_t>& gC) const
{
    if (C.size()!=Z.size()) throw std::invalid_argument("GaussianRF::NuclearGradient: "+std::to_string(C.size())+" nuclei, "+std::to_string(Z.size())+" charges");
    std::vector<double> ea(itsPrims.size()), eb(rb.itsPrims.size());
    for (size_t i=0;i<ea.size();i++) ea[i]=itsPrims[i]->GetExponent();
    for (size_t j=0;j<eb.size();j++) eb[j]=rb.itsPrims[j]->GetExponent();
    const int Ltot=itsL+rb.itsL;
    if (Ω::theNMLs.size()==0) Ω::MakeNMLs();
    const Internal::Boys::AuxillaryFJ Fj(Ltot);   // throws past the Boys table, before the AD call
    const NucPairCtx c{ea.data(), itsCoeff.data(), ea.size(), itsL, eb.data(), rb.itsCoeff.data(), eb.size(), rb.itsL,
                       pa, na, pb, nb, W, Z.data(), Z.size(), &Ω::GetNMLs(Ltot), &Fj};
    std::vector<double> x(6+3*C.size()), g(x.size(),0.0);
    x[0]=itsCenter.x; x[1]=itsCenter.y; x[2]=itsCenter.z; x[3]=rb.itsCenter.x; x[4]=rb.itsCenter.y; x[5]=rb.itsCenter.z;
    for (size_t k=0;k<C.size();k++) { x[6+3*k]=C[k].x; x[7+3*k]=C[k].y; x[8+3*k]=C[k].z; }
    __enzyme_autodiff<void>(reinterpret_cast<void*>(NucPairAt), enzyme_dup, x.data(), g.data(), enzyme_const, &c);
    gA=rvec3_t(g[0],g[1],g[2]); gB=rvec3_t(g[3],g[4],g[5]);
    gC.assign(C.size(), rvec3_t(0,0,0));
    for (size_t k=0;k<C.size();k++) gC[k]=rvec3_t(g[6+3*k],g[7+3*k],g[8+3*k]);
}

} //namespace qchem::BasisSet::Gaussian::Evaluators::PG_Cart_MnD
