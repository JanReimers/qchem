// File: Response/Imp/Probe.C  The amplitude probe and the channel response matrix.
module;
#include <complex>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
module qchem.Response.Probe;
import qchem.Math;    // cos, sin, abs, isfinite, Pi (the Fourier phases)
import qchem.Blaze;

namespace qchem::Response
{

AmplitudeProbe::AmplitudeProbe(const Reference& ref, std::vector<std::vector<cmat_t>> amp, std::vector<std::string> labels)
    : itsRef(ref), itsAmp(std::move(amp)), itsLabels(std::move(labels))
{
    if (itsAmp.size()!=ref.NumBlocks()) throw std::invalid_argument("AmplitudeProbe: one amplitude set per reference block");
    for (size_t b=0;b<itsAmp.size();b++)
    {
        if (itsAmp[b].size()!=itsLabels.size()) throw std::invalid_argument("AmplitudeProbe: one amplitude matrix per channel");
        for (size_t J=0;J<itsLabels.size();J++)
        {
            if (itsAmp[b][J].columns()!=ref.NumOrbitals(b))
                throw std::invalid_argument("AmplitudeProbe: an amplitude matrix whose columns are not the block's orbitals");
            if (itsAmp[b][J].rows()!=itsAmp[0][J].rows())
                throw std::invalid_argument("AmplitudeProbe: a channel with a different function count on different blocks");
        }
    }
}

// A^J on the pair (k, k+q): <psi_{m,k+q}| P_J e^{iq.R} |psi_{n,k}> = sum_mu conj(l_{mu m}(k+q)) l_{mu n}(k).
BlockPairs AmplitudeProbe::Perturbation(size_t J, const MeshShift& q) const
{
    const std::vector<size_t> p=itsRef.Partners(q);
    BlockPairs A;
    A.m.reserve(p.size());
    for (size_t b=0;b<p.size();b++) A.m.push_back(cmat_t(blazem::ctrans(itsAmp[p[b]][J])*itsAmp[b][J]));
    return A;
}

cvec_t AmplitudeProbe::Measure(const MeshShift& q, const BlockPairs& dD) const
{
    cvec_t n(NumChannels());
    for (size_t I=0;I<NumChannels();I++) n[I]=itsRef.Contract(Perturbation(I,q), dD);
    return n;
}

Outcome<ChannelResponse,ResponseFailure> IndependentResponse(const Reference& ref, const ChannelProbe& probe, ivec3_t Nq)
{
    using O=Outcome<ChannelResponse,ResponseFailure>;
    auto qs=ref.QMesh(Nq);
    if (!qs) return O::Fail(qs.Error());
    ChannelResponse r;
    r.Nq=Nq;
    r.gap=std::numeric_limits<double>::infinity();
    r.noise=ref.EigenNoise();
    for (size_t I=0;I<probe.NumChannels();I++) r.labels.push_back(probe.Label(I));
    const size_t nch=probe.NumChannels();
    for (const MeshShift& q : *qs)
    {
        auto g=ref.Gap(q);                          // E1: gate BEFORE any weight is formed
        if (!g) return O::Fail(g.Error());
        r.gap=std::min(r.gap, g->gap);
        cmat_t chi(nch,nch);
        for (size_t J=0;J<nch;J++)
        {
            const cvec_t n=probe.Measure(q, ref.ApplyR0(q, probe.Perturbation(J,q)));
            for (size_t I=0;I<nch;I++) chi(I,J)=n[I];
        }
        r.q.push_back(q);
        r.chi.push_back(chi);
    }
    return O::Ok(std::move(r));
}

rmat_t ChannelResponse::RealSpace() const
{
    const size_t nch=labels.size();
    const size_t NR=size_t(Nq.x)*Nq.y*Nq.z;
    if (q.size()!=NR) throw std::logic_error("ChannelResponse::RealSpace: the q set is not the whole q-mesh");
    std::vector<ivec3_t> R;
    for (int x=0;x<Nq.x;x++) for (int y=0;y<Nq.y;y++) for (int z=0;z<Nq.z;z++) R.push_back(ivec3_t(x,y,z));
    rmat_t out(nch*NR, nch*NR);
    double maxRe=0.0, maxIm=0.0;
    for (size_t a=0;a<NR;a++)
        for (size_t c=0;c<NR;c++)
        {
            const rvec3_t d(R[a].x-R[c].x, R[a].y-R[c].y, R[a].z-R[c].z);
            for (size_t I=0;I<nch;I++)
                for (size_t J=0;J<nch;J++)
                {
                    dcmplx s=0.0;
                    for (size_t iq=0;iq<q.size();iq++)
                    {
                        const rvec3_t qq=q[iq].q();
                        const double ph=2.0*Pi*(qq.x*d.x+qq.y*d.y+qq.z*d.z);
                        s+=dcmplx(cos(ph),sin(ph))*chi[iq](I,J);
                    }
                    s/=double(NR);
                    out(a*nch+I, c*nch+J)=s.real();
                    maxRe=std::max(maxRe,std::abs(s.real()));
                    maxIm=std::max(maxIm,std::abs(s.imag()));
                }
        }
    if (maxIm>1e-10*std::max(maxRe,1e-300))
    {
        std::ostringstream os;
        os << "ChannelResponse::RealSpace: the real-space response is not real (max |Im| " << maxIm << " vs max |Re| "
           << maxRe << ") -- time reversal chi(-q)=conj(chi(q)) is broken";
        throw std::logic_error(os.str());
    }
    return out;
}

std::ostream& ChannelResponse::Write(std::ostream& os, double perUnit, const std::string& unitName) const
{
    const size_t nch=labels.size();
    os << "[chi0] independent-particle channel response, q-mesh " << Nq.x << "x" << Nq.y << "x" << Nq.z
       << ", " << nch << " channels, unit " << unitName << "; response gap "
       << std::setprecision(6) << gap << " Ha (eigenvalue noise " << noise << " Ha";
    if (isfinite(gap) && gap>0) os << ", chi0 relative bound " << std::setprecision(2) << noise/gap;
    os << ")\n";
    for (size_t iq=0;iq<q.size();iq++)
    {
        os << "[chi0]   q=" << q[iq] << "  diag:";
        for (size_t I=0;I<nch;I++) os << " " << std::fixed << std::setprecision(7) << chi[iq](I,I).real()/perUnit;
        os << std::defaultfloat << "\n";
    }
    const rmat_t R=RealSpace();
    os << "[chi0]   real-space, home cell (R=R'=0):\n";
    for (size_t I=0;I<nch;I++)
    {
        os << "[chi0]     " << std::setw(14) << std::left << labels[I] << std::right;
        for (size_t J=0;J<nch;J++) os << " " << std::fixed << std::setprecision(7) << std::setw(11) << R(I,J)/perUnit;
        os << std::defaultfloat << "\n";
    }
    return os;
}

} // namespace
