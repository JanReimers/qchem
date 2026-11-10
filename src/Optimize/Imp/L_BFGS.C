// File: Optimize/Imp/L_BFGS.C  Limited-memory BFGS + trust radius (the header says why).
module;
#include <cmath>
#include <deque>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>
module qchem.Optimize.L_BFGS;
import qchem.Blaze;

namespace qchem::Optimize
{

namespace
{
rvec_t Apply(const Projected* P, const rvec_t& q, const rvec_t& v) {return P ? P->Project(q, v) : v;}
}

L_BFGS::L_BFGS() : L_BFGS(Options{}) {}

L_BFGS::L_BFGS(const Options& options) : itsOpts(options), itsTrust(options.trust)
{
    if (!(itsOpts.initialCurvature>0.0)) throw std::invalid_argument("L_BFGS: the initial curvature must be > 0");
    if (itsOpts.memory==0)               throw std::invalid_argument("L_BFGS: the memory must be > 0");
}

void L_BFGS::Start(const Cost& cost, const rvec_t& q, const CostValue& here)
{
    itsProjected=dynamic_cast<const Projected*>(&cost);
    itsTrust.Reset();
    itsS.clear(); itsY.clear();
    itsQ=q; itsG=here.gradient; itsValue=here.value;
}

rvec_t L_BFGS::InverseHessianTimes(const rvec_t& g) const
{
    const size_t m=itsS.size();
    std::vector<double> alpha(m), rho(m);
    rvec_t r=g;
    for (size_t k=m;k-->0;)
    {
        rho[k]=1.0/blazem::dot(itsY[k], itsS[k]);
        alpha[k]=rho[k]*blazem::dot(itsS[k], r);
        r-=alpha[k]*itsY[k];
    }
    const double gamma = m ? blazem::dot(itsS[m-1], itsY[m-1])/blazem::dot(itsY[m-1], itsY[m-1]) : 1.0/itsOpts.initialCurvature;
    r*=gamma;
    for (size_t k=0;k<m;++k)
    {
        const double beta=rho[k]*blazem::dot(itsY[k], r);
        r+=(alpha[k]-beta)*itsS[k];
    }
    return r;
}

Proposal L_BFGS::Propose(const rvec_t& q, const CostValue& here)
{
    itsQ=q; itsG=here.gradient; itsValue=here.value;
    Proposal p; p.bound=itsTrust.Radius();
    rvec_t d=Apply(itsProjected, q, rvec_t(-InverseHessianTimes(itsG)));
    double gd=blazem::dot(itsG, d);
    if (!(gd<0.0))   // not a descent direction (a pair that broke positive definiteness): forget the memory, steepest descent
    {
        itsS.clear(); itsY.clear();
        d=Apply(itsProjected, q, rvec_t(-itsG/itsOpts.initialCurvature));
        gd=blazem::dot(itsG, d);
    }
    const double l=blazem::norm(d);
    if (!(l>0.0)) return p;   // nothing may move, or nothing to move for: an empty proposal
    double a=1.0;
    if (l>p.bound) { a=p.bound/l; p.limited=true; }
    p.step=a*d;
    // The quadratic model with B = H⁻¹ along d = −H g: m(a d) = a gᵀd + ½ a² dᵀBd = (a − a²/2) gᵀd.
    itsPredicted=(a-0.5*a*a)*gd;
    itsLimited=p.limited;
    return p;
}

Verdict L_BFGS::Tell(const rvec_t& step, const CostValue& there, double noise)
{
    // The pair from the step that was TAKEN, accepted or not (BFGS-RFO's reason: a rejected step's curvature is measured all the same).
    {
        const rvec_t y=Apply(itsProjected, rvec_t(itsQ+step), rvec_t(there.gradient-itsG));
        const double sy=blazem::dot(step, y);
        if (sy>1e-12*blazem::norm(step)*blazem::norm(y))
        {
            itsS.push_back(step); itsY.push_back(y);
            if (itsS.size()>itsOpts.memory) { itsS.pop_front(); itsY.pop_front(); }
        }
    }
    const Verdict v=itsTrust.Judge(there.value-itsValue, itsPredicted, noise, itsLimited);
    if (v.kind!=Verdict::Kind::Accept) return v;
    itsQ+=step; itsG=there.gradient; itsValue=there.value;
    return v;
}

} // namespace qchem::Optimize
