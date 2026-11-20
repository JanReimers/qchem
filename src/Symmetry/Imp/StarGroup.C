// File: Symmetry/Imp/StarGroup.C  The G-space requests of the star-average group (the T1 {G}-star fold).
module;
#include <vector>
#include <complex>
#include <functional>
module qchem.Symmetry.StarGroup;
import qchem.Symmetry.Lattice_3D.Fold;   // FoldGVectors: the star partition of a G set

namespace qchem::Symmetry
{
namespace
{
//! My reciprocal face: \f$\{U|\tau\}\f$ with \f$U=W^\top\f$, the G-index map the fold works in.
std::vector<SymOp> GIndexOps(const std::vector<SymOp>& ops)
{
    std::vector<SymOp> u;
    u.reserve(ops.size());
    for (const SymOp& g : ops) u.push_back({Transpose(g.W), g.tau});
    return u;
}
}

size_t StarGroup::CountGStars(const std::vector<ivec3_t>& gs) const
{
    if (itsOps.empty()) return gs.size();
    return Lattice_3D::FoldGVectors(gs, GIndexOps(itsOps)).Reps();
}

std::vector<dcmplx> StarGroup::EvaluateOnGStars(const std::vector<ivec3_t>& gs, const std::function<dcmplx(const ivec3_t&)>& field,
                                                size_t* nReps) const
{
    std::vector<dcmplx> out(gs.size());
    if (itsOps.empty())
    {
        if (nReps) *nReps=gs.size();
        for (size_t i=0;i<gs.size();i++) out[i]=field(gs[i]);
        return out;
    }
    const std::vector<SymOp> u=GIndexOps(itsOps);
    const Lattice_3D::Fold f=Lattice_3D::FoldGVectors(gs, u);
    if (nReps) *nReps=f.Reps();
    for (size_t r=0;r<f.Reps();++r)
    {
        const ivec3_t& mr=gs[f.repRaw[r]];
        const dcmplx   v =field(mr);
        for (auto [mi, o] : f.members[r])
            if (o<0) out[mi]=field(gs[mi]);                       // no edge op: evaluate directly (still exact)
            else
            {
                const double ph=2.0*Pi*(mr.x*u[o].tau.x + mr.y*u[o].tau.y + mr.z*u[o].tau.z);
                out[mi]=std::polar(1.0, ph)*v;                    // f(U m_rep) = e^{+2πi m_rep·τ} f(m_rep)
            }
    }
    return out;
}

} // namespace qchem::Symmetry
