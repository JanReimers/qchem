// File: GaussianRF_AD.C  THE CARTESIAN GAUSSIAN DERIVATIVES, Enzyme's (GeometryOptimizationPlan, "DERIVATIVES ARE ENZYME'S"): the basis
// gradient in r, and the 2-, 3- and 4-centre kernels' geometry derivatives (KetGradient, NuclearGradient, ThreeCentreGradient,
// FourCentreGradient, below).
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
#include <memory>
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

// The 2-centre COULOMB (a|b) (TwoCentre::Repulsion: the Coulomb fit's metric), the cache-free twin of PrimGaussian::Repulsion2C: per
// primitive pair the Gaussian product, F_j(alpha |A-B|^2) with alpha = ab/p and the Hermite-Coulomb recursion at A-B on a flat store (what
// Ω::SelfRNLM caches), contracted by the shared PairRepulsion2C.  The Hermite1 coefficients of both sides depend on the exponents alone, so
// they are tabulated OUTSIDE the differentiated function (h: per bra primitive i and component p at hOff[i*na+p], per ket primitive j and
// component q at hOff[nA*na+j*nb+q], each a box over NLM <= the component) and read as constants.
struct KetCoulombCtx
{
    const double* ea; const double* ca; size_t nA; int LA; const double* A;
    const double* eb; const double* cb; size_t nB; int LB;
    const Polarization* pa; size_t na; const Polarization* pb; size_t nb;
    const double* W;
    const std::vector<Polarization>* NLMs; const Internal::Boys::AuxillaryFJ* Fj;
    const double* h; const size_t* hOff;
};
double KetCoulombAt(const double* x, const KetCoulombCtx* c)
{
    using Internal::MnD::Triangle3D;
    const rvec3_t A(c->A[0],c->A[1],c->A[2]), B(x[0],x[1],x[2]);
    const int Lr=c->LA+c->LB;   // as Ω::Ltotal
    auto tri=[](int n){return size_t(n+1)*(n+2)*(n+3)/6;};
    std::vector<size_t> off(Lr+2,0);
    for (int j=0;j<=Lr;j++) off[j+1]=off[j]+tri(Lr-j);
    std::vector<double> store(off[Lr+1],0.0), Fj(Lr+1,0.0);
    auto at=[&](int j, int N, int L, int M) -> double&
    {   const int n=Lr-j; return store[off[j]+Triangle3D::Layer(n-N)+Triangle3D::Row(n-N-L)+size_t(n-N-L-M)]; };
    auto Rof=[&](const Polarization& m){ return store[Triangle3D::Layer(Lr-m.n)+Triangle3D::Row(Lr-m.n-m.l)+size_t(Lr-m.n-m.l-m.m)]; };
    auto box=[](const double* t, const Polarization& P){ return [t,&P](const Polarization& n){ return t[(size_t(n.n)*(P.l+1)+n.l)*(P.m+1)+n.m]; }; };
    double s=0.0;
    for (size_t i=0;i<c->nA;i++)
        for (size_t j=0;j<c->nB;j++)
        {
            const GaussProduct gp(c->ea[i],A,c->eb[j],B);
            const double alpha=gp.ab/gp.p;   // as Ω::SelfRNLM
            c->Fj->GetFjAt(alpha*(gp.AB*gp.AB), Fj);
            Internal::MnD::HermiteCoulomb(Lr, alpha, gp.AB, Fj, at);
            double t=0.0;
            for (size_t p=0;p<c->na;p++)
                for (size_t q=0;q<c->nb;q++)
                {
                    const double w=c->W[p*c->nb+q];
                    if (w==0.0) continue;   // a constant weight
                    const double* ha=c->h+c->hOff[i*c->na+p];
                    const double* hb=c->h+c->hOff[c->nA*c->na+j*c->nb+q];
                    t += w*PairRepulsion2C(gp.ab,gp.p,*c->NLMs,c->pa[p],c->pb[q],box(ha,c->pa[p]),box(hb,c->pb[q]),Rof);
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
    double x[3]={B.x,B.y,B.z}, g[3]={0.0,0.0,0.0};
    if (K==TwoCentre::Repulsion)
    {
        std::vector<double> h; std::vector<size_t> hOff;
        auto table=[&](const GaussianRF& r, const Polarization* P, size_t nP)
        {
            for (const auto& prim : r.itsPrims)
            {
                const Hermite1& H1=prim->GetH1();
                for (size_t u=0;u<nP;u++)
                {
                    hOff.push_back(h.size());
                    for (int n=0;n<=P[u].n;n++) for (int l=0;l<=P[u].l;l++) for (int m=0;m<=P[u].m;m++)
                        h.push_back(H1(Polarization(n,l,m),P[u]));
                }
            }
        };
        table(*this, pa, na);
        table(rb, pb, nb);
        const Internal::Boys::AuxillaryFJ Fj(itsL+rb.itsL);   // throws past the Boys table, before the AD call
        if (Ω::theNMLs.size()==0) Ω::MakeNMLs();
        const KetCoulombCtx c{ea.data(), itsCoeff.data(), ea.size(), itsL, A, eb.data(), rb.itsCoeff.data(), eb.size(), rb.itsL,
                              pa, na, pb, nb, W, &Ω::GetNMLs(itsL), &Fj, h.data(), hOff.data()};
        __enzyme_autodiff<void>(reinterpret_cast<void*>(KetCoulombAt), enzyme_dup, x, g, enzyme_const, &c);
        return rvec3_t(g[0],g[1],g[2]);
    }
    const KetPairCtx c{K, ea.data(), itsCoeff.data(), ea.size(), itsL, A, eb.data(), rb.itsCoeff.data(), eb.size(), rb.itsL,
                       pa, na, pb, nb, W};
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


// THE 3-CENTRE KERNELS' GEOMETRY DERIVATIVE, every centre moving (lane B 2026-11-13; it replaced libcint's int3c1e_ip / int3c2e_ip1/_ip2
// for the molecular PP field pullbacks).  The cache-free twins of PrimGaussian::Overlap3C (a Hermite3 built from the shared
// MakeThreeCentreProduct per primitive triple) and PrimGaussian::Repulsion3C (the pair's Hermite2, the triple's Hermite-Coulomb table on a
// flat store, the shared PairRepulsion3C).  c's Hermite1 coefficients depend on its exponent alone, so they are tabulated OUTSIDE the
// differentiated function (hc: per primitive k, term t, a box over NLM' <= p_t) and read as constants.  x = [A, B, C].
struct ThreeCCtx
{
    GaussianRF::ThreeCentre K;
    const double* ea; const double* ca; size_t nA; int LA;
    const double* eb; const double* cb; size_t nB; int LB;
    const double* ec; const double* cc; size_t nCp; int LC;
    const Polarization* pa; size_t na; const Polarization* pb; size_t nb;
    const Polarization* pt; const double* ct; size_t nt;
    const double* W;
    const std::vector<Polarization>* NLMs; const Internal::Boys::AuxillaryFJ* Fj;
    const double* hc; const size_t* hcOff;   // hc[hcOff[k*nt+t] + (n*(p_t.l+1)+l)*(p_t.m+1)+m]
};
double ThreeCAt(const double* x, const ThreeCCtx* c)
{
    using Internal::MnD::Triangle3D;
    const rvec3_t A(x[0],x[1],x[2]), B(x[3],x[4],x[5]), C(x[6],x[7],x[8]);
    double s=0.0;
    if (c->K==GaussianRF::ThreeCentre::Overlap)
    {
        for (size_t i=0;i<c->nA;i++)
            for (size_t j=0;j<c->nB;j++)
                for (size_t k=0;k<c->nCp;k++)
                {
                    const ThreeCentreProduct tp=MakeThreeCentreProduct(c->ea[i],A,c->eb[j],B,c->ec[k],C);
                    const Hermite3 H3(tp.alphaQ,tp.QA,tp.QB,tp.QC,c->LA,c->LB,c->LC,tp.Eabc);
                    double t=0.0;
                    for (size_t p=0;p<c->na;p++)
                        for (size_t q=0;q<c->nb;q++)
                        {
                            const double w=c->W[p*c->nb+q];
                            if (w==0.0) continue;   // a constant weight
                            for (size_t u=0;u<c->nt;u++) t += w*c->ct[u]*H3(c->pa[p],c->pb[q],c->pt[u]);
                        }
                    s += c->ca[i]*c->cb[j]*c->cc[k]*t;
                }
        return s;
    }
    const int Lab=c->LA+c->LB, Lr=Lab+c->LC;
    auto tri=[](int n){return size_t(n+1)*(n+2)*(n+3)/6;};
    std::vector<size_t> off(Lr+2,0);
    for (int j=0;j<=Lr;j++) off[j+1]=off[j]+tri(Lr-j);
    std::vector<double> store(off[Lr+1],0.0), Fj(Lr+1,0.0);
    auto at=[&](int j, int N, int L, int M) -> double&
    {   const int n=Lr-j; return store[off[j]+Triangle3D::Layer(n-N)+Triangle3D::Row(n-N-L)+size_t(n-N-L-M)]; };
    auto Rof=[&](const Polarization& m){ return store[Triangle3D::Layer(Lr-m.n)+Triangle3D::Row(Lr-m.n-m.l)+size_t(Lr-m.n-m.l-m.m)]; };
    for (size_t i=0;i<c->nA;i++)
        for (size_t j=0;j<c->nB;j++)
        {
            const GaussProduct gp(c->ea[i],A,c->eb[j],B);
            const Hermite2 H2(gp.p,gp.PA,gp.PB,c->LA+1,c->LB+1);   // as Ω's
            for (size_t k=0;k<c->nCp;k++)
            {
                const double alpha=gp.p*c->ec[k]/(gp.p+c->ec[k]);   // as findRNLM
                const rvec3_t dR=gp.P-C;
                c->Fj->GetFjAt(alpha*dR*dR, Fj);
                Internal::MnD::HermiteCoulomb(Lr, alpha, dR, Fj, at);
                double t=0.0;
                for (size_t p=0;p<c->na;p++)
                    for (size_t q=0;q<c->nb;q++)
                    {
                        const double w=c->W[p*c->nb+q];
                        if (w==0.0) continue;
                        for (size_t u=0;u<c->nt;u++)
                        {
                            const Polarization& pu=c->pt[u];
                            const double* h=c->hc+c->hcOff[k*c->nt+u];
                            auto Hc=[&](const Polarization& q3){ return h[(size_t(q3.n)*(pu.l+1)+q3.l)*(pu.m+1)+q3.m]; };
                            t += w*c->ct[u]*PairRepulsion3C(H2,gp.p,gp.Eij,c->ec[k],c->pa[p],c->pb[q],pu,*c->NLMs,Hc,Rof);
                        }
                    }
                s += c->ca[i]*c->cb[j]*c->cc[k]*t;
            }
        }
    return s;
}

void GaussianRF::ThreeCentreGradient(ThreeCentre K, rf_t& ra, rf_t& rb, const Polarization* pa, size_t na, const Polarization* pb,
                                     size_t nb, const Polarization* pt, const double* ct, size_t nt, const double* W,
                                     rvec3_t& gA, rvec3_t& gB, rvec3_t& gC) const
{
    auto exps=[](const GaussianRF& r){ std::vector<double> e(r.itsPrims.size()); for (size_t i=0;i<e.size();i++) e[i]=r.itsPrims[i]->GetExponent(); return e; };
    const std::vector<double> ea=exps(ra), eb=exps(rb), ec=exps(*this);
    std::vector<double> hc; std::vector<size_t> hcOff;
    std::unique_ptr<Internal::Boys::AuxillaryFJ> Fj;
    if (K==ThreeCentre::Repulsion)
    {
        for (size_t k=0;k<itsPrims.size();k++)
        {
            const Hermite1& H1=itsPrims[k]->GetH1();
            for (size_t u=0;u<nt;u++)
            {
                hcOff.push_back(hc.size());
                for (int n=0;n<=pt[u].n;n++) for (int l=0;l<=pt[u].l;l++) for (int m=0;m<=pt[u].m;m++)
                    hc.push_back(H1(Polarization(n,l,m),pt[u]));
            }
        }
        Fj=std::make_unique<Internal::Boys::AuxillaryFJ>(ra.itsL+rb.itsL+itsL);   // throws past the Boys table, before the AD call
        if (Ω::theNMLs.size()==0) Ω::MakeNMLs();
    }
    else if (Ω::theNMLs.size()==0) Ω::MakeNMLs();
    const ThreeCCtx c{K, ea.data(), ra.itsCoeff.data(), ea.size(), ra.itsL, eb.data(), rb.itsCoeff.data(), eb.size(), rb.itsL,
                      ec.data(), itsCoeff.data(), ec.size(), itsL, pa, na, pb, nb, pt, ct, nt, W,
                      &Ω::GetNMLs(ra.itsL+rb.itsL), Fj.get(), hc.data(), hcOff.data()};
    double x[9]={ra.itsCenter.x,ra.itsCenter.y,ra.itsCenter.z, rb.itsCenter.x,rb.itsCenter.y,rb.itsCenter.z, itsCenter.x,itsCenter.y,itsCenter.z};
    double g[9]={0,0,0,0,0,0,0,0,0};
    __enzyme_autodiff<void>(reinterpret_cast<void*>(ThreeCAt), enzyme_dup, x, g, enzyme_const, &c);
    gA=rvec3_t(g[0],g[1],g[2]); gB=rvec3_t(g[3],g[4],g[5]); gC=rvec3_t(g[6],g[7],g[8]);
}


// THE 4-CENTRE COULOMB'S GEOMETRY DERIVATIVE, every centre moving (lane B 2026-11-14; it replaced libcint's int2e_ip1 for the molecular HF
// J and K pullbacks).  The cache-free twin of PrimGaussian::Repulsion4C: per primitive pair the Gaussian product and its Hermite2 (as Ω
// builds them), per primitive quartet F_j(alpha |P-Q|^2) with alpha = pq/(p+q) and the Hermite-Coulomb recursion at P-Q on a flat store
// (what findRNLM caches), contracted by the shared PairRepulsion4C with the shared QuartetLambda.  x = [A, B, C, D]: one reverse pass
// gives all four centres.  The cd pair's product is rebuilt per ab pair (no table of Hermite2 objects inside the differentiated
// function).  Gated by tests/M_KetGradient.C (the FD of the cached kernel) and tests/M_GeometryDerivatives.C (the FD of the J / K
// contractions; libcint's ERI4Gradient as an oracle).
struct FourCCtx
{
    const double* e[4]; const double* c[4]; size_t np[4]; int L[4];
    const Polarization* p[4]; size_t n[4];
    const double* W;
    const Internal::Boys::AuxillaryFJ* Fj;
};
double FourCAt(const double* x, const FourCCtx* c)
{
    using Internal::MnD::Triangle3D;
    const rvec3_t A(x[0],x[1],x[2]), B(x[3],x[4],x[5]), C(x[6],x[7],x[8]), D(x[9],x[10],x[11]);
    const int Lr=c->L[0]+c->L[1]+c->L[2]+c->L[3];   // as findRNLM: Ltotal(ab) + Ltotal(cd)
    auto tri=[](int n){return size_t(n+1)*(n+2)*(n+3)/6;};
    std::vector<size_t> off(Lr+2,0);
    for (int j=0;j<=Lr;j++) off[j+1]=off[j]+tri(Lr-j);
    std::vector<double> store(off[Lr+1],0.0), Fj(Lr+1,0.0);
    auto at=[&](int j, int N, int L, int M) -> double&
    {   const int n=Lr-j; return store[off[j]+Triangle3D::Layer(n-N)+Triangle3D::Row(n-N-L)+size_t(n-N-L-M)]; };
    auto R=[&](int N, int L, int M){ return store[Triangle3D::Layer(Lr-N)+Triangle3D::Row(Lr-N-L)+size_t(Lr-N-L-M)]; };
    const size_t na=c->n[0], nb=c->n[1], nc=c->n[2], nd=c->n[3];
    double s=0.0;
    for (size_t i=0;i<c->np[0];i++)
        for (size_t j=0;j<c->np[1];j++)
        {
            const GaussProduct ab(c->e[0][i],A,c->e[1][j],B);
            const Hermite2 Hab(ab.p,ab.PA,ab.PB,c->L[0]+1,c->L[1]+1);   // as Ω's
            for (size_t k=0;k<c->np[2];k++)
                for (size_t l=0;l<c->np[3];l++)
                {
                    const GaussProduct cd(c->e[2][k],C,c->e[3][l],D);
                    const Hermite2 Hcd(cd.p,cd.PA,cd.PB,c->L[2]+1,c->L[3]+1);
                    const double alpha=ab.p*cd.p/(ab.p+cd.p);   // as findRNLM
                    const rvec3_t dR=ab.P-cd.P;
                    c->Fj->GetFjAt(alpha*(dR*dR), Fj);
                    Internal::MnD::HermiteCoulomb(Lr, alpha, dR, Fj, at);
                    const double lambda=QuartetLambda(ab.p,ab.Eij,cd.p,cd.Eij);
                    double t=0.0;
                    for (size_t a=0;a<na;a++)
                        for (size_t b=0;b<nb;b++)
                            for (size_t q=0;q<nc;q++)
                                for (size_t d=0;d<nd;d++)
                                {
                                    const double w=c->W[((a*nb+b)*nc+q)*nd+d];
                                    if (w==0.0) continue;   // a constant weight
                                    t += w*PairRepulsion4C(Hab,c->p[0][a],c->p[1][b],Hcd,c->p[2][q],c->p[3][d],lambda,R);
                                }
                    s += c->c[0][i]*c->c[1][j]*c->c[2][k]*c->c[3][l]*t;
                }
        }
    return s;
}

void GaussianRF::FourCentreGradient(rf_t& ra, rf_t& rb, rf_t& rc, const Polarization* pa, size_t na, const Polarization* pb, size_t nb,
                                    const Polarization* pc, size_t nc, const Polarization* pd, size_t nd, const double* W,
                                    rvec3_t& gA, rvec3_t& gB, rvec3_t& gC, rvec3_t& gD) const
{
    auto exps=[](const GaussianRF& r){ std::vector<double> e(r.itsPrims.size()); for (size_t i=0;i<e.size();i++) e[i]=r.itsPrims[i]->GetExponent(); return e; };
    const std::vector<double> ea=exps(ra), eb=exps(rb), ec=exps(rc), ed=exps(*this);
    const Internal::Boys::AuxillaryFJ Fj(ra.itsL+rb.itsL+rc.itsL+itsL);   // throws past the Boys table, before the AD call
    const FourCCtx c{{ea.data(), eb.data(), ec.data(), ed.data()},
                     {ra.itsCoeff.data(), rb.itsCoeff.data(), rc.itsCoeff.data(), itsCoeff.data()},
                     {ea.size(), eb.size(), ec.size(), ed.size()}, {ra.itsL, rb.itsL, rc.itsL, itsL},
                     {pa, pb, pc, pd}, {na, nb, nc, nd}, W, &Fj};
    double x[12]={ra.itsCenter.x,ra.itsCenter.y,ra.itsCenter.z, rb.itsCenter.x,rb.itsCenter.y,rb.itsCenter.z,
                  rc.itsCenter.x,rc.itsCenter.y,rc.itsCenter.z, itsCenter.x,itsCenter.y,itsCenter.z};
    double g[12]={0,0,0,0,0,0,0,0,0,0,0,0};
    __enzyme_autodiff<void>(reinterpret_cast<void*>(FourCAt), enzyme_dup, x, g, enzyme_const, &c);
    gA=rvec3_t(g[0],g[1],g[2]); gB=rvec3_t(g[3],g[4],g[5]); gC=rvec3_t(g[6],g[7],g[8]); gD=rvec3_t(g[9],g[10],g[11]);
}

} //namespace qchem::BasisSet::Gaussian::Evaluators::PG_Cart_MnD
