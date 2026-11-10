// File: Optimize/Imp/BFGS_RFO.C  BFGS + RFO + trust radius (the header says why, and the three noise rules).
module;
#include <cmath>
#include <cstdio>
#include <optional>
#include <stdexcept>
#include <string>
module qchem.Optimize.BFGS_RFO;
import qchem.Blaze;

namespace qchem::Optimize
{

namespace
{
std::string Sci(double x) {char b[32]; std::snprintf(b, sizeof b, "%.3g", x); return b;}
rvec_t Apply(const Projected* P, const rvec_t& q, const rvec_t& v) {return P ? P->Project(q, v) : v;}

//! An orthonormal basis (the columns) of the step space at q: the eigenvectors of the projector with eigenvalue ~ 1.
rmat_t StepBasis(const Projected* P, const rvec_t& q)
{
    const size_t n=q.size();
    if (!P) { rmat_t B(n, n, 0.0); for (size_t i=0;i<n;++i) B(i,i)=1.0; return B; }
    rhmat_t M(n);
    for (size_t j=0;j<n;++j)
    {
        rvec_t e(n, 0.0); e[j]=1.0;
        const rvec_t c=P->Project(q, e);
        for (size_t i=0;i<=j;++i) M(i,j)=c[i];   // P = Pᵀ (an orthogonal projector): the upper triangle carries it
    }
    rvec_t w; rmat_t U;
    blazem::eigen(M, w, U);
    size_t m=0; for (double x : w) if (x>0.5) ++m;
    rmat_t B(n, m);
    size_t k=0;
    for (size_t j=0;j<n;++j) if (w[j]>0.5) { for (size_t i=0;i<n;++i) B(i,k)=U(i,j); ++k; }
    return B;
}
}

BFGS_RFO::BFGS_RFO() : BFGS_RFO(Options{}) {}

BFGS_RFO::BFGS_RFO(const Options& options) : itsOpts(options)
{
    if (!(itsOpts.trustInitial>0.0 && itsOpts.trustMin>0.0 && itsOpts.trustMax>=itsOpts.trustInitial))
        throw std::invalid_argument("BFGS_RFO: the trust radii must satisfy 0 < trustMin, 0 < trustInitial <= trustMax");
    if (!(itsOpts.initialCurvature>0.0)) throw std::invalid_argument("BFGS_RFO: the initial curvature must be > 0");
}

void BFGS_RFO::Start(const Cost& cost, const rvec_t& q, const CostValue& here)
{
    const size_t n=q.size();
    itsProjected=dynamic_cast<const Projected*>(&cost);
    itsH=rmat_t(n, n, 0.0);
    if (const auto* M=dynamic_cast<const ModelHessian*>(&cost))
    {
        for (size_t j=0;j<n;++j) { rvec_t e(n, 0.0); e[j]=1.0; const rvec_t c=M->ApplyH0(q, e); for (size_t i=0;i<n;++i) itsH(i,j)=c[i]; }
        itsH=0.5*(itsH+blazem::trans(itsH));
    }
    else
        for (size_t i=0;i<n;++i) itsH(i,i)=itsOpts.initialCurvature;
    itsTrust=itsOpts.trustInitial;
    itsQ=q; itsG=here.gradient; itsValue=here.value;
}

Proposal BFGS_RFO::Propose(const rvec_t& q, const CostValue& here)
{
    itsQ=q; itsG=here.gradient; itsValue=here.value;
    // THE RFO STEP in the step space: the lowest eigenvector of [[Hr, gr], [grᵀ, 0]] gives s = v/v_last, a descent direction for any Hr.
    const rmat_t B=StepBasis(itsProjected, q);
    const size_t m=B.columns();
    Proposal p; p.bound=itsTrust;
    if (m==0) return p;   // nothing may move
    const rmat_t Hr=blazem::trans(B)*itsH*B;
    const rvec_t gr=blazem::trans(B)*itsG;
    rhmat_t A(m+1);
    for (size_t i=0;i<m;++i) { for (size_t j=i;j<m;++j) A(i,j)=0.5*(Hr(i,j)+Hr(j,i)); A(i,m)=gr[i]; }
    A(m,m)=0.0;
    rvec_t w; rmat_t V;
    blazem::eigen(A, w, V);
    rvec_t sr(m);
    if (std::fabs(V(m,0))>1e-12) for (size_t i=0;i<m;++i) sr[i]=V(i,0)/V(m,0);
    else                         sr=-gr;   // the RFO's degenerate case: fall back on steepest descent, which the radius then caps
    if (const double l=blazem::norm(sr); l>itsTrust) { sr*=itsTrust/l; p.limited=true; }
    p.step=Apply(itsProjected, q, rvec_t(B*sr));
    itsPredicted=blazem::dot(itsG, p.step)+0.5*blazem::dot(p.step, rvec_t(itsH*p.step));
    itsLimited=p.limited;
    return p;
}

Verdict BFGS_RFO::Tell(const rvec_t& step, const CostValue& there, double noise)
{
    // BFGS on the step that was TAKEN, accepted or not: a rejected step's (s, y) is curvature measured all the same, and dropping it
    // wastes the evaluation (the 6-d bowl from 1.5 bohr: 3 rejections without this update, 2 with it, 2026-11-07).
    // H += y yᵀ/(yᵀs) − (Hs)(Hs)ᵀ/(sᵀHs), skipped when the curvature condition fails, which keeps H positive definite.
    {
        const rvec_t q1=itsQ+step;
        const rvec_t y=Apply(itsProjected, q1, rvec_t(there.gradient-itsG));
        const double sy=blazem::dot(step, y);
        if (sy>1e-12*blazem::norm(step)*blazem::norm(y))
        {
            const rvec_t Hs=itsH*step;
            const double sHs=blazem::dot(step, Hs);
            itsH+=(y*blazem::trans(y))/sy;
            if (sHs>0.0) itsH-=(Hs*blazem::trans(Hs))/sHs;
        }
    }
    const double change=there.value-itsValue;
    if (change>noise)   // RULE 1: only a rise beyond the noise is a rise.  Undo the step, cut the radius.
    {
        itsTrust*=0.25;
        if (itsTrust<itsOpts.trustMin)
            return {Verdict::Kind::Collapsed, "the trust radius fell below "+Sci(itsOpts.trustMin)+" after a value rise of "
                                              +Sci(change)+" (noise bound "+Sci(noise)+")"};
        return {Verdict::Kind::Reject, {}};
    }
    // RULE 2: the ratio means something only when the model predicted a change above the noise; below it, hold the radius.
    if (std::fabs(itsPredicted)>noise)
    {
        const double ratio = itsPredicted<0.0 ? change/itsPredicted : 0.0;
        if      (ratio>0.75 && itsLimited) itsTrust=std::min(2.0*itsTrust, itsOpts.trustMax);
        else if (ratio<0.25)               itsTrust=std::max(0.5*itsTrust, itsOpts.trustMin);
    }
    itsQ+=step; itsG=there.gradient; itsValue=there.value;
    return {Verdict::Kind::Accept, {}};
}

} // namespace qchem::Optimize
