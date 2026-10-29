// File: Response/QStars.C  ONE q PER STAR (doc/KAxisPlan.md P3; doc/Records/LinearResponsePlan.md Q9).
//
// The q-mesh response computes chi_q for every q, but symmetry-equivalent q are related exactly.  With ToRealSpace's
// convention chi(I,R_a; J,R_c) = (1/N) Sum_q e^{2 pi i q.(R_a-R_c)} X_q(I,J), and an op g = {W|tau} of the GROUND STATE's
// group that maps site I to site gI shifted by the lattice vector t_I (W x_I + tau = x_{gI} + t_I, fractional), the
// invariance chi(g(I,R_a); g(J,R_c)) = chi(I,R_a; J,R_c) gives
//
//     X_{q'}(gI, gJ) = e^{-2 pi i q'.(t_I - t_J)} X_q(I, J),      q' = W^{-T} q  (mod G).
//
// So one solved q per star gives every member by a channel PERMUTATION and a PHASE -- no orbital rotation, for channels
// that are scalar per manifold (HubbardResolution::Manifold, hp.x's).  A spin-flip op (sigma) leaves a spin-summed
// manifold channel unchanged, so magnetic (Shubnikov) ops serve too.  Using a group ASSERTS the ground state has it
// (KAxisPlan §5 ruling 1): the caller opts in, and the response AUDITS one star member by solving it directly.
module;
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>
export module qchem.Response.QStars;
export import qchem.Symmetry.Lattice_3D.BlochQN;   // MeshShift
export import qchem.Symmetry.SymOp;               // SymOp {W, tau, sigma}
import qchem.Matrix3D;
import qchem.Types;

export namespace qchem::Response
{
using Symmetry::Lattice_3D::MeshShift;

//! What the group acts on, per response channel: the SITE it lives on and its angular momentum (a channel maps to the
//! channel of the same l on the image site).
struct ChannelSite { size_t site=0; int l=0; };

//! How q-mesh point i is obtained: from star representative \c rep, under the op that sends channel c to \c perm[c] with
//! the site's lattice shift \c shift[c].  A representative is its own member (identity permutation, zero shifts).
struct StarMember
{
    size_t              rep=0;
    std::vector<size_t> perm;
    std::vector<ivec3_t> shift;
};

//! The q-mesh folded into stars (one \c StarMember per q, in \c Reference::QMesh order).
struct QStars
{
    std::vector<StarMember> member;
    bool   IsRep  (size_t i) const {return member[i].rep==i;}
    size_t NumReps()         const {size_t n=0; for (size_t i=0;i<member.size();i++) n+=IsRep(i); return n;}
};

//! \brief Fold the q-mesh \a q under \a ops (the GROUND STATE's group, fractional direct frame).  \a siteFrac: every site's
//! fractional position; \a channels: every response channel's (site, l); \a perturbed: the perturbed channels, which an op
//! must map onto themselves to be used.  An op that sends a channel to a (site, l) with no channel is skipped, as is one
//! that does not map the mesh onto itself.  THROWS on a site an op cannot place (a group that is not the structure's).
QStars FoldQMesh(const std::vector<MeshShift>& q, const std::vector<Symmetry::SymOp>& ops,
                 const std::vector<rvec3_t>& siteFrac, const std::vector<ChannelSite>& channels,
                 const std::vector<size_t>& perturbed, double tol=1e-6)
{
    const size_t nq=q.size(), nc=channels.size();
    auto frac=[](double x){return x-std::floor(x+1e-9);};                       // into [0,1), robust at 1-eps
    auto apply=[](const Matrix3D<double>& M, const rvec3_t& v)                   // Matrix3D is 1-based
    { return rvec3_t(M(1,1)*v.x+M(1,2)*v.y+M(1,3)*v.z, M(2,1)*v.x+M(2,2)*v.y+M(2,3)*v.z, M(3,1)*v.x+M(3,2)*v.y+M(3,3)*v.z); };
    QStars out;
    out.member.assign(nq, StarMember{});
    std::vector<bool> placed(nq,false);
    for (size_t i=0;i<nq;i++)
    {
        if (placed[i]) continue;
        StarMember self; self.rep=i; self.perm.resize(nc); self.shift.assign(nc, ivec3_t(0,0,0));
        for (size_t c=0;c<nc;c++) self.perm[c]=c;
        out.member[i]=self; placed[i]=true;
        for (const auto& g : ops)
        {
            // The channel action: site s -> s' with W x_s + tau = x_{s'} + t.
            std::vector<size_t> perm(nc); std::vector<ivec3_t> shift(nc);
            bool usable=true;
            for (size_t c=0;c<nc && usable;c++)
            {
                const rvec3_t y=apply(g.W, siteFrac[channels[c].site])+g.tau;
                size_t img=siteFrac.size();
                for (size_t s=0;s<siteFrac.size();s++)
                {
                    const rvec3_t d=y-siteFrac[s];
                    if (std::abs(d.x-std::round(d.x))<tol && std::abs(d.y-std::round(d.y))<tol && std::abs(d.z-std::round(d.z))<tol)
                    { img=s; shift[c]=ivec3_t(int(std::round(d.x)),int(std::round(d.y)),int(std::round(d.z))); break; }
                }
                if (img==siteFrac.size())
                    throw std::invalid_argument("FoldQMesh: an op does not map a site onto a site -- not the structure's group");
                size_t cImg=nc;
                for (size_t k=0;k<nc;k++) if (channels[k].site==img && channels[k].l==channels[c].l) {cImg=k; break;}
                if (cImg==nc) usable=false; else perm[c]=cImg;
            }
            for (size_t J : perturbed)                                           // the perturbed set must map onto itself
            {
                bool in=false; for (size_t K : perturbed) in = in || perm[J]==K;
                usable = usable && in;
            }
            if (!usable) continue;
            // q' = W^{-T} q, reduced onto the mesh.
            const rvec3_t qp=apply(Transpose(Invert(g.W)), q[i].q());
            const rvec3_t r(frac(qp.x), frac(qp.y), frac(qp.z));
            for (size_t j=0;j<nq;j++)
            {
                const rvec3_t d=r-q[j].q();
                if (std::abs(d.x-std::round(d.x))<tol && std::abs(d.y-std::round(d.y))<tol && std::abs(d.z-std::round(d.z))<tol)
                {
                    if (!placed[j]) { out.member[j]={i, perm, shift}; placed[j]=true; }
                    break;
                }
            }
        }
    }
    return out;
}

} //namespace
