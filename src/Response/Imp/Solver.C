// File: Response/Imp/Solver.C  The self-consistent linear response by GMRES.
module;
#include <optional>
#include <algorithm>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>
module qchem.Response.Solver;
import qchem.Math;    // Pi (the q-star phases)
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
        const BlockPairs R0KdD=itsRef.ApplyR0(*itsRule, InducedFockMO(itsFrame, itsK, dD, itsRule));
        return x-itsRef.Pack(R0KdD);
    }
private:
    const Reference&                                 itsRef;
    const OrbitalFrame<T>&                           itsFrame;
    const Hamiltonian::ResponseKernel<T>&            itsK;
    std::shared_ptr<const Symmetry::SelectionRule>   itsRule;
    size_t                                           itsN;
};

//! The rule a little-group reduction hands the kernel: \a inner's pairs whose KET is an orbit representative (by its spatial
//! block identity -- both spins of a representative k are gathered), with \a inner's wave vector.
class RepresentativeKets : public virtual Symmetry::SelectionRule, public virtual Symmetry::WaveVectorShift
{
public:
    RepresentativeKets(std::shared_ptr<const Symmetry::SelectionRule> inner, std::vector<size_t> kets)
        : itsInner(std::move(inner)), itsKets(std::move(kets)) {}
    virtual bool Couples(const Symmetry::Symmetry& bra, const Symmetry::Symmetry& ket) const override
    {
        return itsInner->Couples(bra, ket)
            && std::find(itsKets.begin(), itsKets.end(), ket.SequenceIndex())!=itsKets.end();
    }
    virtual rvec3_t q() const override {return ChargeDensity::WaveVectorOf(*itsInner);}
private:
    std::shared_ptr<const Symmetry::SelectionRule> itsInner;
    std::vector<size_t>                            itsKets;
};

//! \f$x\mapsto x-\mathcal R_0\mathcal K\,{\rm Herm}(x)\f$ on a LITTLE-GROUP packing (doc/KAxisPlan.md P5): x is unfolded onto every
//! block (the scatter builds the whole δρ), the kernel gathers δF on the orbit representatives only, and δF is unfolded
//! again for R0 (whose q = 0 Fermi shift sums over every block).  The kernel and its terms never see the group.
class ReducedResponseOperator : public LinearOperator<dcmplx>
{
public:
    ReducedResponseOperator(const Reference& ref, const OrbitalFrame<dcmplx>& frame, const Hamiltonian::ResponseKernel<dcmplx>& K,
                            std::shared_ptr<const Symmetry::SelectionRule> rule, const LittleGroupReduction& red, size_t n)
        : itsRef(ref), itsFrame(frame), itsK(K), itsRule(std::move(rule)), itsRed(red), itsN(n)
    {
        std::vector<size_t> kets;
        for (size_t b : red.Representatives()) kets.push_back(ref.BlockIrrep(b).sym->SequenceIndex());
        itsReps=std::make_shared<RepresentativeKets>(itsRule, std::move(kets));
    }
    virtual size_t Dimension() const override {return itsN;}
    virtual cvec_t Apply(const cvec_t& x, double) const override
    {
        BlockPairs dD=itsRef.HermitianPart(itsRed.Unpack(x), *itsRule);
        double s=0.0;                                   // InducedFockMO's unit-scale guard (see its declaration)
        for (const auto& m : dD.m)
            for (size_t i=0;i<m.rows();i++) for (size_t j=0;j<m.columns();j++) s=std::max(s, std::abs(m(i,j)));
        if (s==0.0) return x;
        for (auto& m : dD.m) m/=s;
        const auto delta=ChargeDensity::RestrictCoupling(itsFrame.ToAO(dD, itsRule), itsReps);
        const auto dF=itsK.InducedFock(*delta);
        BlockPairs X=itsRed.Unfold(itsFrame.ToMO(*dF, *itsRule, itsRed.Representatives()));
        for (auto& m : X.m) m*=s;
        return x-itsRed.Pack(itsRef.ApplyR0(*itsRule, X));
    }
private:
    const Reference&                                 itsRef;
    const OrbitalFrame<dcmplx>&                      itsFrame;
    const Hamiltonian::ResponseKernel<dcmplx>&       itsK;
    std::shared_ptr<const Symmetry::SelectionRule>   itsRule;
    std::shared_ptr<const Symmetry::SelectionRule>   itsReps;
    const LittleGroupReduction&                      itsRed;
    size_t                                           itsN;
};

} // namespace

template <class T> BlockPairs InducedFockMO(const OrbitalFrame<T>& frame, const Hamiltonian::ResponseKernel<T>& K, const BlockPairs& dD,
                                            std::shared_ptr<const Symmetry::SelectionRule> rule)
{
    double s=0.0;
    for (const auto& m : dD.m)
        for (size_t i=0;i<m.rows();i++) for (size_t j=0;j<m.columns();j++) s=std::max(s, std::abs(m(i,j)));
    BlockPairs X;
    if (s==0.0)                                   // K[0] = 0: the zero vector, shaped as ToMO would give it
    {
        for (const auto& m : dD.m) X.m.push_back(cmat_t(m.rows(), m.columns(), dcmplx(0.0)));
        return X;
    }
    BlockPairs unit=dD;
    for (auto& m : unit.m) m/=s;
    const auto delta=frame.ToAO(unit, rule);
    const auto dF=K.InducedFock(*delta);
    X=frame.ToMO(*dF, *rule);
    for (auto& m : X.m) m*=s;
    return X;
}
template BlockPairs InducedFockMO<double>(const OrbitalFrame<double>&, const Hamiltonian::ResponseKernel<double>&,
                                          const BlockPairs&, std::shared_ptr<const Symmetry::SelectionRule>);
template BlockPairs InducedFockMO<dcmplx>(const OrbitalFrame<dcmplx>&, const Hamiltonian::ResponseKernel<dcmplx>&,
                                          const BlockPairs&, std::shared_ptr<const Symmetry::SelectionRule>);

template <class T> Outcome<SelfConsistentResponse,ResponseFailure>
LinearResponse(const Reference& ref, const OrbitalFrame<T>& frame, const Hamiltonian::ResponseKernel<T>& kernel,
               const ChannelProbe& probe, std::shared_ptr<const Symmetry::SelectionRule> rule, const KrylovParams& kp,
               std::vector<size_t> perturbed, const OrbitalImages* images)
{
    if constexpr (!std::is_same_v<T,dcmplx>)
        if (images) throw std::logic_error("LinearResponse: a little-group reduction needs complex (Bloch) blocks");
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
        // THE KRYLOV TABLE (user, 2026-09-29: a 47-minute stall that printed nothing is unreadable): one row per
        // operator application -- the recurrence's residual estimate -- and the TRUE residual at every restart.
        KrylovParams kpc=kp;
        if (!kpc.observer)
        {
            std::cout << "[krylov] channel " << probe.Label(perturbed[c]) << ": GMRES(" << kp.restart << "), tol "
                      << kp.tol << ", budget " << kp.maxIter << " kernel applications" << std::endl
                      << "[krylov]      #   residual" << std::endl;
            kpc.observer=[](const KrylovStep& st)
            {
                std::cout << "[krylov] " << std::setw(6) << st.iteration << "   " << std::scientific << std::setprecision(3)
                          << st.residual << std::defaultfloat << (st.trueResidual ? "   (TRUE, restart)" : "") << std::endl;
            };
        }
        // THE LITTLE GROUP of (q, J) (KAxisPlan P5): solve on its orbit representatives, measure on the unfolded δD.
        std::optional<LittleGroupReduction> red;
        if constexpr (std::is_same_v<T,dcmplx>)
            if (images)
            {
                red.emplace(ref, *images, *rule, probe.Perturbation(perturbed[c], *rule));
                std::cout << "[little group] channel " << probe.Label(perturbed[c]) << ": ";
                red->Write(std::cout);
            }
        auto solve=[&]()
        {
            if constexpr (std::is_same_v<T,dcmplx>)
                if (red)
                {
                    const cvec_t br=red->Pack(R0V[c]);
                    ReducedResponseOperator Ar(ref, frame, kernel, rule, *red, br.size());
                    return SolveGMRES<dcmplx>(Ar, br, nullptr, kpc);
                }
            return SolveGMRES<dcmplx>(A, b, nullptr, kpc);
        };
        auto s=solve();
        if (!s)
        {
            std::ostringstream os;
            os << "channel " << probe.Label(perturbed[c]) << ": " << s.Error().detail;
            return O::Fail({ResponseFailure::Why::NotConverged, os.str()});
        }
        r.residual  .push_back(s->residual);
        r.iterations.push_back(s->iterations);
        const cvec_t n0=probe.Measure(*rule, R0V[c]);
        const cvec_t n1=probe.Measure(*rule, red ? red->Unpack(s->x) : ref.Unpack(s->x, *rule));
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
        // 8 significant figures: enough to hold a chi against the finite-difference cross-check (the stream's own
        // precision may have been left at 3 by an earlier writer -- it printed "-2.07" against FD -2.0672).
        const std::streamsize prec=os.precision(8);
        os << "[response]   " << std::setw(8) << labels[J] << "  chi0 " << std::setw(14) << chi0(J,c).real()
           << "  chi " << std::setw(14) << chi(J,c).real() << "  residual " << std::setprecision(3) << residual[c]
           << " (" << iterations[c] << " kernel applications)" << std::endl;
        os.precision(prec);
    }
    if (perturbed.size()>1)       // the J x J block, off-diagonals included: what a supercell U inverts, so a LOG can rebuild U
    {
        const std::streamsize prec=os.precision(10);
        for (const char* which : {"chi0","chi"})
        {
            const cmat_t M = std::string(which)=="chi0" ? Chi0JJ() : ChiJJ();
            for (size_t a=0;a<M.rows();a++)
            {
                os << "[response]   " << std::setw(4) << which << " JxJ row " << std::setw(8) << labels[perturbed[a]] << ":";
                for (size_t c=0;c<M.columns();c++) os << "  (" << M(a,c).real() << "," << M(a,c).imag() << ")";
                os << std::endl;
            }
        }
        os.precision(prec);
    }
    return os;
}

template <class T> Outcome<QMeshResponse,ResponseFailure>
LinearResponse(const Reference& ref, const OrbitalFrame<T>& frame, const Hamiltonian::ResponseKernel<T>& kernel,
               const ChannelProbe& probe, ivec3_t Nq, const KrylovParams& kp, std::vector<size_t> perturbed,
               const QStars* stars, double starTol, const OrbitalImages* images)
{
    using O=Outcome<QMeshResponse,ResponseFailure>;
    auto qs=ref.QMesh(Nq);
    if (!qs) return O::Fail(qs.Error());
    if (stars && stars->member.size()!=qs->size()) throw std::logic_error("LinearResponse: the q-star fold is not over this q-mesh");
    QMeshResponse r;
    r.Nq=Nq;
    const size_t nq=qs->size();
    std::vector<std::optional<SelfConsistentResponse>> solved(nq);
    // THE AUDIT member (KAxisPlan §5 ruling 1): the first non-representative q, solved directly as well.
    size_t audit=nq;
    if (stars) for (size_t i=0;i<nq && audit==nq;i++) if (!stars->IsRep(i)) audit=i;
    if (stars) std::cout << "[response] q-stars: " << stars->NumReps() << " of " << nq << " q solved (one per star)"
                         << (audit<nq ? ", plus one member as the symmetry audit" : "") << std::endl;
    for (size_t i=0;i<nq;i++)
    {
        const MeshShift& q=(*qs)[i];
        if (stars && !stars->IsRep(i) && i!=audit) continue;
        std::cout << "[response] q = " << q << "  (" << i+1 << " of " << nq << ")" << (i==audit ? "  [audit]" : "") << std::endl;
        auto s=LinearResponse(ref, frame, kernel, probe, std::make_shared<MeshShift>(q), kp, perturbed, images);
        if (!s)
        {
            std::ostringstream os;
            os << "q = " << q << ": " << s.Error().detail;
            return O::Fail({s.Error().why, os.str()});
        }
        s->Write(std::cout);                     // reported at its own activity (pin 17), before the next q starts
        solved[i]=std::move(*s);
    }
    // GENERATE the other members: X_{q'}(perm I, perm J) = e^{-2 pi i q'.(t_I - t_J)} X_q(I, J)  (QStars.C).
    auto generate=[&](size_t i)
    {
        const StarMember& m=stars->member[i];
        const SelfConsistentResponse& src=*solved[m.rep];
        SelfConsistentResponse out=src;
        const rvec3_t qq=(*qs)[i].q();
        const std::vector<size_t>& P=src.perturbed;
        auto col=[&](size_t channel){ for (size_t k=0;k<P.size();k++) if (P[k]==channel) return k;
                                      throw std::logic_error("LinearResponse: a star op leaves the perturbed set"); };
        for (size_t I=0;I<src.chi0.rows();I++)
            for (size_t J=0;J<P.size();J++)
            {
                const ivec3_t dt=m.shift[I]-m.shift[P[J]];
                const double ph=-2.0*Pi*(qq.x*dt.x+qq.y*dt.y+qq.z*dt.z);
                const dcmplx e(std::cos(ph), std::sin(ph));
                out.chi0(m.perm[I], col(m.perm[P[J]]))=e*src.chi0(I,J);
                out.chi (m.perm[I], col(m.perm[P[J]]))=e*src.chi (I,J);
            }
        return out;
    };
    if (audit<nq)
    {
        const SelfConsistentResponse g=generate(audit), &d=*solved[audit];
        double worst=0.0, scale=0.0;
        for (size_t I=0;I<d.chi0.rows();I++) for (size_t J=0;J<d.chi0.columns();J++)
        {
            worst=std::max({worst, std::abs(g.chi0(I,J)-d.chi0(I,J)), std::abs(g.chi(I,J)-d.chi(I,J))});
            scale=std::max({scale, std::abs(d.chi0(I,J)), std::abs(d.chi(I,J))});
        }
        std::cout << "[response] q-star audit at q = " << (*qs)[audit] << ": max|generated - solved| / max|chi| = "
                  << worst/scale << "  (tolerance " << starTol << ")" << std::endl;
        if (!(worst<=starTol*scale))
        {
            std::ostringstream os;
            os << "the q-star audit FAILED at q = " << (*qs)[audit] << " (relative " << worst/scale << " > " << starTol
               << "): the ground state does not have the group the star fold assumed -- run without qStars";
            return O::Fail({ResponseFailure::Why::Configuration, os.str()});
        }
    }
    // Generate EVERY member before moving any representative out (a generated member reads its representative).
    std::vector<SelfConsistentResponse> all(nq);
    for (size_t i=0;i<nq;i++) if (!solved[i]) all[i]=generate(i);
    for (size_t i=0;i<nq;i++)
    {
        r.q.push_back((*qs)[i]);
        r.atQ.push_back(solved[i] ? std::move(*solved[i]) : std::move(all[i]));
    }
    return O::Ok(std::move(r));
}

namespace {
//! The time-reversal tolerance of a SOLVED chi's real-space transform: the worst Krylov residual, times the conditioning
//! headroom of (1 - R0 K) (measured ~190x on NiO k222), floored at chi0's own 1e-8.
double ChiRealnessTol(const QMeshResponse& r)
{
    double res=0.0;
    for (const auto& s : r.atQ) for (double x : s.residual) res=std::max(res, x);
    return std::max(1e-8, 1e4*res);
}
} // namespace

rmat_t QMeshResponse::Chi0RealSpace(double* imagRatio) const
{
    std::vector<cmat_t> X;
    for (const auto& s : atQ) X.push_back(s.Chi0JJ());
    return ToRealSpace(Nq, q, X, 1e-8, imagRatio);
}
rmat_t QMeshResponse::ChiRealSpace(double* imagRatio) const
{
    std::vector<cmat_t> X;
    for (const auto& s : atQ) X.push_back(s.ChiJJ());
    return ToRealSpace(Nq, q, X, ChiRealnessTol(*this), imagRatio);
}

rmat_t QMeshResponse::Interaction() const
{
    return rmat_t(blazem::inv(Chi0RealSpace())-blazem::inv(ChiRealSpace()));
}

rvec_t QMeshResponse::UScalar() const
{
    const rmat_t X0=Chi0RealSpace(), X=ChiRealSpace();
    const size_t nJ=atQ.empty() ? 0 : atQ[0].perturbed.size();
    rvec_t u(nJ);
    for (size_t c=0;c<nJ;c++) u[c]=1.0/X0(c,c)-1.0/X(c,c);
    return u;
}

rvec_t QMeshResponse::U() const
{
    const rmat_t W=Interaction();
    const size_t nJ=atQ.empty() ? 0 : atQ[0].perturbed.size();
    rvec_t u(nJ);
    for (size_t c=0;c<nJ;c++) u[c]=W(c,c);      // the home cell, R = 0: indices 0 .. nJ-1
    return u;
}

std::ostream& QMeshResponse::Write(std::ostream& os) const
{
    if (atQ.empty()) return os;
    const auto& lab=atQ[0].labels;
    const auto& J  =atQ[0].perturbed;
    const size_t nJ=J.size();
    double leak0=0.0, leak=0.0;
    const rmat_t X0=Chi0RealSpace(&leak0), X=ChiRealSpace(&leak);
    const std::streamsize prec=os.precision(8);
    os << "[response] real space over the " << Nq.x << "x" << Nq.y << "x" << Nq.z << " q-mesh's supercell, home cell "
       << "(R = R' = 0), perturbed x perturbed, 1/Ha:" << std::endl;
    for (size_t a=0;a<nJ;a++)
    {
        os << "[response]   " << std::setw(14) << std::left << lab[J[a]] << std::right << "  chi0";
        for (size_t c=0;c<nJ;c++) os << " " << std::setw(14) << X0(a,c);
        os << "   chi";
        for (size_t c=0;c<nJ;c++) os << " " << std::setw(14) << X(a,c);
        os << std::endl;
    }
    os << "[response]   time-reversal leak max|Im|/max|Re|: chi0 " << std::setprecision(2) << leak0 << ", chi " << leak
       << std::setprecision(8) << std::endl;
    const rvec_t u=U();
    os << "[response] U = (chi0^-1 - chi^-1) over the " << Nq.x << "x" << Nq.y << "x" << Nq.z << " q-mesh's supercell "
       << "(+U frozen; " << nJ << " perturbed channels):";
    for (size_t c=0;c<nJ;c++) os << "  " << lab[J[c]] << " " << std::setprecision(6) << u[c]*27.211386245988 << " eV";
    os << std::endl;
    if (nJ>1)                                                   // with one channel every inversion is the same number
    {
        const rvec_t us=UScalar();
        os << "[response] U 1x1 = 1/chi0_cc - 1/chi_cc (each channel alone in the supercell):";
        for (size_t c=0;c<nJ;c++) os << "  " << lab[J[c]] << " " << std::setprecision(6) << us[c]*27.211386245988 << " eV";
        os << std::endl;
    }
    os.precision(prec);
    return os;
}

LR_HubbardU::LR_HubbardU(QMeshResponse r, std::vector<Hamiltonian::HubbardChannel> channels)
    : itsR(std::move(r)), itsChannels(std::move(channels))
{
    if (itsR.atQ.empty()) throw std::invalid_argument("LR_HubbardU: a response with no q point");
    if (itsR.atQ[0].labels.size()!=itsChannels.size())
        throw std::invalid_argument("LR_HubbardU: the response's channels are not the term's manifolds");
}

std::vector<Hamiltonian::HubbardEstimate> LR_HubbardU::Evaluate() const
{
    const rvec_t U=itsR.U();
    const auto& J=itsR.atQ[0].perturbed;
    std::vector<Hamiltonian::HubbardEstimate> out;
    for (size_t c=0;c<J.size();c++) out.push_back({.manifold=J[c], .site=itsChannels[J[c]].site, .l=itsChannels[J[c]].l, .U=U[c]});
    return out;
}

std::ostream& LR_HubbardU::Write(std::ostream& os) const
{
    const std::streamsize prec=os.precision(6);
    for (const auto& e : Evaluate())
        os << "[LR] site " << e.site << " l=" << e.l << ": U=" << e.U*27.211386245988 << " eV  ((chi0^-1 - chi^-1) over the "
           << itsR.Nq.x << "x" << itsR.Nq.y << "x" << itsR.Nq.z << " q-mesh's supercell, +U frozen)" << std::endl;
    os.precision(prec);
    return os;
}

template Outcome<SelfConsistentResponse,ResponseFailure> LinearResponse<double>(const Reference&, const OrbitalFrame<double>&,
    const Hamiltonian::ResponseKernel<double>&, const ChannelProbe&, std::shared_ptr<const Symmetry::SelectionRule>, const KrylovParams&, std::vector<size_t>, const OrbitalImages*);
template Outcome<SelfConsistentResponse,ResponseFailure> LinearResponse<dcmplx>(const Reference&, const OrbitalFrame<dcmplx>&,
    const Hamiltonian::ResponseKernel<dcmplx>&, const ChannelProbe&, std::shared_ptr<const Symmetry::SelectionRule>, const KrylovParams&, std::vector<size_t>, const OrbitalImages*);

template Outcome<QMeshResponse,ResponseFailure> LinearResponse<double>(const Reference&, const OrbitalFrame<double>&,
    const Hamiltonian::ResponseKernel<double>&, const ChannelProbe&, ivec3_t, const KrylovParams&, std::vector<size_t>, const QStars*, double, const OrbitalImages*);
template Outcome<QMeshResponse,ResponseFailure> LinearResponse<dcmplx>(const Reference&, const OrbitalFrame<dcmplx>&,
    const Hamiltonian::ResponseKernel<dcmplx>&, const ChannelProbe&, ivec3_t, const KrylovParams&, std::vector<size_t>, const QStars*, double, const OrbitalImages*);

} // namespace
