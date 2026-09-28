// File: Response/Imp/Reference.C  The independent-particle response R0 by sum over states.
module;
#include <complex>   // std::conj on a dcmplx scalar (NOT blaze::conj -- a no-op on scalars)
#include <iomanip>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
module qchem.Response.Reference;
import qchem.Streamable;   // Irrep's operator<< (the failure messages name the block)
import qchem.Math;         // fabs, isfinite
import qchem.Blaze;

namespace qchem::Response
{

Reference::Reference(std::vector<ReferenceBlock> blocks, std::unique_ptr<OccupancyRule> rule, double eigenNoise)
    : itsBlocks(std::move(blocks)), itsRule(std::move(rule)), itsNoise(eigenNoise)
{
    if (!itsRule) throw std::invalid_argument("Response::Reference: no occupancy rule");
    if (itsBlocks.empty()) throw std::invalid_argument("Response::Reference: no blocks");
    for (const auto& b : itsBlocks)
    {
        if (b.e.size()!=b.f.size()) throw std::invalid_argument("Response::Reference: eigenvalue/occupation length mismatch");
        // D5: one stored block per mesh point.  A star > 1 is an IBZ representative standing for points
        // that are NOT stored, and a k+q partner may be one of them.
        if (b.irrep.sym->GetDegeneracy()!=1)
            throw std::invalid_argument("Response::Reference: an IBZ-reduced block (k-star > 1) -- a linear response "
                "needs the FULL k-mesh stored (doc/LinearResponsePlan.md D5); run without symmetry imposition");
    }
}

Outcome<std::vector<MeshShift>,ResponseFailure> Reference::QMesh(ivec3_t Nq) const
{
    using O=Outcome<std::vector<MeshShift>,ResponseFailure>;
    auto qs=Symmetry::Lattice_3D::CommensurateShifts(*itsBlocks[0].irrep.sym, Nq);
    if (!qs) return O::Fail({ResponseFailure::Why::Incommensurate, qs.Error()});
    return O::Ok(qs.TakeValue());
}

std::vector<size_t> Reference::Partners(const SelectionRule& rule) const
{
    std::vector<size_t> out(itsBlocks.size());
    for (size_t b=0;b<itsBlocks.size();b++)
    {
        size_t n=0;
        for (size_t c=0;c<itsBlocks.size();c++)
            if (itsBlocks[c].irrep.ms==itsBlocks[b].irrep.ms &&
                rule.Couples(*itsBlocks[c].irrep.sym, *itsBlocks[b].irrep.sym)) {out[b]=c; n++;}
        if (n!=1)
        {
            std::ostringstream os;
            os << "Response::Reference::Partners: block " << itsBlocks[b].irrep << " has " << n
               << " stored partners under this selection rule (exactly 1 needed; for a lattice, a full, unreduced"
                  " k-mesh, D5)";
            throw std::logic_error(os.str());
        }
    }
    return out;
}

Outcome<ResponseGap,ResponseFailure> Reference::Gap(const SelectionRule& rule) const
{
    using O=Outcome<ResponseGap,ResponseFailure>;
    ResponseGap r;
    r.noise=itsNoise;
    if (!itsRule->RequiresResolvedGap()) return O::Ok(r);   // a smeared weight is bounded: nothing to gate
    const std::vector<size_t> p=Partners(rule);
    size_t bw=0, nw=0, mw=0;
    for (size_t b=0;b<itsBlocks.size();b++)
    {
        const ReferenceBlock& k=itsBlocks[b];
        const ReferenceBlock& kq=itsBlocks[p[b]];
        for (size_t n=0;n<k.e.size();n++)
            for (size_t m=0;m<kq.e.size();m++)
            {
                if (k.f[n]==kq.f[m]) continue;
                // Δ = energy of the LESS occupied level minus that of the MORE occupied one: > 0 is aufbau order.
                const double d = k.f[n]>kq.f[m] ? kq.e[m]-k.e[n] : k.e[n]-kq.e[m];
                if (d<r.gap) {r.gap=d; bw=b; nw=n; mw=m;}
            }
    }
    if (r.gap<=0.0 || r.gap<=itsNoise)   // a NaN noise compares false: only the sign gates
    {
        const ReferenceBlock& k=itsBlocks[bw];
        const ReferenceBlock& kq=itsBlocks[p[bw]];
        std::ostringstream os;
        os << std::setprecision(6)
           << (r.gap<=0.0 ? "INVERTED coupled pair" : "UNRESOLVED coupled pair")
           << ": ket block " << k.irrep << " orbital " << nw << " (e=" << k.e[nw] << ", f=" << k.f[nw] << ")"
           << " -> bra block " << kq.irrep << " orbital " << mw << " (e=" << kq.e[mw] << ", f=" << kq.f[mw] << ")"
           << ", gap " << r.gap << " Ha against eigenvalue noise ";
        if (isfinite(itsNoise)) os << itsNoise << " Ha.";
        else                    os << "UNMEASURED (this recipe computes no [F,D]; only the gap's sign was gated).";
        os
           << (r.gap<=0.0 ? "  The per-block integer fill is not an aufbau state ACROSS this pair: the state is"
                            " not a gapped insulator on this mesh -- treat it as a metal (Fermi occupancy)."
                          : "  The gap is not resolved above the reference's own eigenvalue noise.");
        return O::Fail({r.gap<=0.0 ? ResponseFailure::Why::Inverted : ResponseFailure::Why::Unresolved, os.str()});
    }
    return O::Ok(r);
}

BlockPairs Reference::ApplyR0(const SelectionRule& rule, const BlockPairs& dF) const
{
    if (dF.m.size()!=itsBlocks.size()) throw std::invalid_argument("Response::Reference::ApplyR0: one matrix per block");
    const std::vector<size_t> p=Partners(rule);
    // The Fermi shift exists only when the perturbation has a DIAGONAL (δF_nn on one block): every block is
    // its own partner -- q = 0, or a totally symmetric molecular perturbation.  A k -> k+q or A1 -> B1
    // perturbation moves no level at first order, so it cannot move μ.
    bool selfPaired=true;
    for (size_t b=0;b<p.size();b++) selfPaired = selfPaired && p[b]==b;
    BlockPairs dD;
    dD.m.resize(itsBlocks.size());
    // The q = 0 Fermi shift, per reservoir: the numerator and denominator of δμ_r.
    std::vector<dcmplx> num;
    std::vector<double> den;
    auto grow=[&](int r){ if (r>=int(num.size())) {num.resize(r+1,0.0); den.resize(r+1,0.0);} };
    for (size_t b=0;b<itsBlocks.size();b++)
    {
        const ReferenceBlock& k =itsBlocks[b];
        const ReferenceBlock& kq=itsBlocks[p[b]];
        if (k.g!=kq.g) throw std::logic_error("Response::Reference::ApplyR0: a pair of blocks with different level capacities");
        const cmat_t& F=dF.m[b];
        if (F.rows()!=kq.e.size() || F.columns()!=k.e.size())
            throw std::invalid_argument("Response::Reference::ApplyR0: a block pair's matrix has the wrong shape");
        cmat_t& X=dD.m[b];
        X.resize(F.rows(), F.columns());
        for (size_t m=0;m<kq.e.size();m++)
            for (size_t n=0;n<k.e.size();n++)
                X(m,n) = k.g*itsRule->ResponseWeight(k.e[n],k.f[n],kq.e[m],kq.f[m]) * F(m,n);
        if (selfPaired)
        {
            grow(k.reservoir);
            for (size_t n=0;n<k.e.size();n++)
            {
                const double D=k.g*itsRule->ResponseWeight(k.e[n],k.f[n],k.e[n],k.f[n]);   // g f'_n (0 for integer)
                num[k.reservoir]+=k.w*D*F(n,n);
                den[k.reservoir]+=k.w*D;
            }
        }
    }
    if (selfPaired)
        for (size_t b=0;b<itsBlocks.size();b++)
        {
            const ReferenceBlock& k=itsBlocks[b];
            if (den[k.reservoir]==0.0) continue;   // an integer reservoir, or no level at the Fermi level
            const dcmplx dmu=num[k.reservoir]/den[k.reservoir];
            for (size_t n=0;n<k.e.size();n++)
                dD.m[b](n,n) -= k.g*itsRule->ResponseWeight(k.e[n],k.f[n],k.e[n],k.f[n])*dmu;
        }
    return dD;
}

dcmplx Reference::Contract(const BlockPairs& a, const BlockPairs& x) const
{
    if (a.m.size()!=itsBlocks.size() || x.m.size()!=itsBlocks.size())
        throw std::invalid_argument("Response::Reference::Contract: one matrix per block");
    dcmplx s=0.0;
    for (size_t b=0;b<itsBlocks.size();b++)
    {
        const cmat_t& A=a.m[b];
        const cmat_t& X=x.m[b];
        dcmplx t=0.0;
        for (size_t i=0;i<A.rows();i++)
            for (size_t j=0;j<A.columns();j++) t+=std::conj(A(i,j))*X(i,j);
        s+=itsBlocks[b].w*t;
    }
    return s;
}

} // namespace
