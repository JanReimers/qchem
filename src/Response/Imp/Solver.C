// File: Response/Imp/Solver.C  The self-consistent linear response by GMRES.
module;
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>
module qchem.Response.Solver;
import qchem.Blaze;

namespace qchem::Response
{

namespace {

//! \f$x\mapsto x-\mathcal R_0\mathcal K\,{\rm Herm}(x)\f$ on the packed orbital-basis δD (the Reference's packing, Q3).
//! The kernel sees only the HERMITIAN part of x (Reference::HermitianPart): the anti-Hermitian part is decoupled
//! -- the operator is the identity on it and the right-hand side has none of it -- so the solution is exactly
//! the one without the projection, while rounding amplified by Gram-Schmidt in late Krylov vectors can no
//! longer reach the kernel as a spurious non-Hermitian δD.
template <class T> class ResponseOperator : public LinearOperator<dcmplx>
{
public:
    ResponseOperator(const Reference& ref, const OrbitalFrame<T>& frame, const Hamiltonian::ResponseKernel<T>& K,
                     std::shared_ptr<const Symmetry::SelectionRule> rule, size_t n)
        : itsRef(ref), itsFrame(frame), itsK(K), itsRule(std::move(rule)), itsN(n) {}
    virtual size_t Dimension() const override {return itsN;}
    virtual cvec_t Apply(const cvec_t& x, double) const override   // exact: the tolerance is not needed
    {
        const BlockPairs dD=itsRef.HermitianPart(itsRef.Unpack(x, *itsRule), *itsRule);
        const auto delta=itsFrame.ToAO(dD, itsRule);
        const auto dF=itsK.InducedFock(*delta);
        const BlockPairs R0KdD=itsRef.ApplyR0(*itsRule, itsFrame.ToMO(*dF, *itsRule));
        return x-itsRef.Pack(R0KdD);
    }
private:
    const Reference&                                 itsRef;
    const OrbitalFrame<T>&                           itsFrame;
    const Hamiltonian::ResponseKernel<T>&            itsK;
    std::shared_ptr<const Symmetry::SelectionRule>   itsRule;
    size_t                                           itsN;
};

} // namespace

template <class T> Outcome<SelfConsistentResponse,ResponseFailure>
LinearResponse(const Reference& ref, const OrbitalFrame<T>& frame, const Hamiltonian::ResponseKernel<T>& kernel,
               const ChannelProbe& probe, std::shared_ptr<const Symmetry::SelectionRule> rule, const KrylovParams& kp)
{
    using O=Outcome<SelfConsistentResponse,ResponseFailure>;
    auto g=ref.Gap(*rule);                          // E1: gate BEFORE any weight is formed
    if (!g) return O::Fail(g.Error());
    SelfConsistentResponse r;
    r.gap=g->gap;
    r.noise=g->noise;
    const size_t nch=probe.NumChannels();
    for (size_t I=0;I<nch;I++) r.labels.push_back(probe.Label(I));
    r.chi0.resize(nch,nch);
    r.chi .resize(nch,nch);

    if (nch==0) return O::Ok(std::move(r));
    std::vector<BlockPairs> R0V(nch);
    for (size_t J=0;J<nch;J++) R0V[J]=ref.ApplyR0(*rule, probe.Perturbation(J, *rule));
    const size_t n=ref.Pack(R0V[0]).size();
    ResponseOperator<T> A(ref, frame, kernel, rule, n);
    for (size_t J=0;J<nch;J++)
    {
        const cvec_t b=ref.Pack(R0V[J]);
        auto s=SolveGMRES<dcmplx>(A, b, nullptr, kp);
        if (!s)
        {
            std::ostringstream os;
            os << "channel " << probe.Label(J) << ": " << s.Error().detail;
            return O::Fail({ResponseFailure::Why::NotConverged, os.str()});
        }
        r.residual  .push_back(s->residual);
        r.iterations.push_back(s->iterations);
        const cvec_t n0=probe.Measure(*rule, R0V[J]);
        const cvec_t n1=probe.Measure(*rule, ref.Unpack(s->x, *rule));
        for (size_t I=0;I<nch;I++) {r.chi0(I,J)=n0[I]; r.chi(I,J)=n1[I];}
    }
    return O::Ok(std::move(r));
}

std::ostream& SelfConsistentResponse::Write(std::ostream& os) const
{
    os << "[response] " << labels.size() << " channels, response gap " << gap << " Ha";
    if (std::isfinite(noise)) os << " (eigenvalue noise " << noise << " Ha)";
    else                      os << " (eigenvalue noise UNMEASURED: only the gap's sign was gated)";
    os << std::endl;
    for (size_t I=0;I<labels.size();I++)
    {
        os << "[response]   " << std::setw(8) << labels[I] << "  chi0 " << std::setw(14) << chi0(I,I).real()
           << "  chi " << std::setw(14) << chi(I,I).real() << "  residual " << residual[I]
           << " (" << iterations[I] << " kernel applications)" << std::endl;
    }
    return os;
}

template Outcome<SelfConsistentResponse,ResponseFailure> LinearResponse<double>(const Reference&, const OrbitalFrame<double>&,
    const Hamiltonian::ResponseKernel<double>&, const ChannelProbe&, std::shared_ptr<const Symmetry::SelectionRule>, const KrylovParams&);
template Outcome<SelfConsistentResponse,ResponseFailure> LinearResponse<dcmplx>(const Reference&, const OrbitalFrame<dcmplx>&,
    const Hamiltonian::ResponseKernel<dcmplx>&, const ChannelProbe&, std::shared_ptr<const Symmetry::SelectionRule>, const KrylovParams&);

} // namespace
