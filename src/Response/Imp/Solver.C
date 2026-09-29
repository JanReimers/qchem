// File: Response/Imp/Solver.C  The self-consistent linear response by GMRES.
module;
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
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
               const ChannelProbe& probe, std::shared_ptr<const Symmetry::SelectionRule> rule, const KrylovParams& kp,
               std::vector<size_t> perturbed)
{
    using O=Outcome<SelfConsistentResponse,ResponseFailure>;
    auto g=ref.Gap(*rule);                          // E1: gate BEFORE any weight is formed
    if (!g) return O::Fail(g.Error());
    SelfConsistentResponse r;
    r.gap=g->gap;
    r.noise=g->noise;
    const size_t nch=probe.NumChannels();
    for (size_t I=0;I<nch;I++) r.labels.push_back(probe.Label(I));
    if (perturbed.empty()) for (size_t J=0;J<nch;J++) perturbed.push_back(J);
    for (size_t J : perturbed)
        if (J>=nch) throw std::out_of_range("LinearResponse: a perturbed channel index beyond the probe's channels");
    r.perturbed=perturbed;
    const size_t nJ=perturbed.size();
    r.chi0.resize(nch,nJ);
    r.chi .resize(nch,nJ);

    if (nch==0) return O::Ok(std::move(r));
    std::vector<BlockPairs> R0V(nJ);
    for (size_t c=0;c<nJ;c++) R0V[c]=ref.ApplyR0(*rule, probe.Perturbation(perturbed[c], *rule));
    const size_t n=ref.Pack(R0V[0]).size();
    ResponseOperator<T> A(ref, frame, kernel, rule, n);
    for (size_t c=0;c<nJ;c++)
    {
        const cvec_t b=ref.Pack(R0V[c]);
        auto s=SolveGMRES<dcmplx>(A, b, nullptr, kp);
        if (!s)
        {
            std::ostringstream os;
            os << "channel " << probe.Label(perturbed[c]) << ": " << s.Error().detail;
            return O::Fail({ResponseFailure::Why::NotConverged, os.str()});
        }
        r.residual  .push_back(s->residual);
        r.iterations.push_back(s->iterations);
        const cvec_t n0=probe.Measure(*rule, R0V[c]);
        const cvec_t n1=probe.Measure(*rule, ref.Unpack(s->x, *rule));
        for (size_t I=0;I<nch;I++) {r.chi0(I,c)=n0[I]; r.chi(I,c)=n1[I];}
    }
    return O::Ok(std::move(r));
}

namespace {
cmat_t RowsOf(const cmat_t& X, const std::vector<size_t>& rows)
{
    cmat_t Y(rows.size(), X.columns());
    for (size_t a=0;a<rows.size();a++) for (size_t c=0;c<X.columns();c++) Y(a,c)=X(rows[a],c);
    return Y;
}
} // namespace
cmat_t SelfConsistentResponse::Chi0JJ() const {return RowsOf(chi0, perturbed);}
cmat_t SelfConsistentResponse::ChiJJ () const {return RowsOf(chi , perturbed);}

std::ostream& SelfConsistentResponse::Write(std::ostream& os) const
{
    os << "[response] " << labels.size() << " channels, response gap " << gap << " Ha";
    if (std::isfinite(noise)) os << " (eigenvalue noise " << noise << " Ha)";
    else                      os << " (eigenvalue noise UNMEASURED: only the gap's sign was gated)";
    os << std::endl;
    if (perturbed.size()<labels.size())
        os << "[response]   perturbed " << perturbed.size() << " of " << labels.size() << " channels (every channel is measured)" << std::endl;
    for (size_t c=0;c<perturbed.size();c++)
    {
        const size_t J=perturbed[c];
        os << "[response]   " << std::setw(8) << labels[J] << "  chi0 " << std::setw(14) << chi0(J,c).real()
           << "  chi " << std::setw(14) << chi(J,c).real() << "  residual " << residual[c]
           << " (" << iterations[c] << " kernel applications)" << std::endl;
    }
    return os;
}

template Outcome<SelfConsistentResponse,ResponseFailure> LinearResponse<double>(const Reference&, const OrbitalFrame<double>&,
    const Hamiltonian::ResponseKernel<double>&, const ChannelProbe&, std::shared_ptr<const Symmetry::SelectionRule>, const KrylovParams&, std::vector<size_t>);
template Outcome<SelfConsistentResponse,ResponseFailure> LinearResponse<dcmplx>(const Reference&, const OrbitalFrame<dcmplx>&,
    const Hamiltonian::ResponseKernel<dcmplx>&, const ChannelProbe&, std::shared_ptr<const Symmetry::SelectionRule>, const KrylovParams&, std::vector<size_t>);

} // namespace
