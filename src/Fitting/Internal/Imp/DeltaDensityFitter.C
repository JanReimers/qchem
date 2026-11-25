// File: Fitting/Internal/Imp/DeltaDensityFitter.C  The bodies of qchem.Fitting.Internal.DeltaDensityFitter: a density fit on a δ fit basis, its
// Hartree energy and its fixed-D force legs (the interface carries the design and the rulings).
module;
#include <functional>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
module qchem.Fitting.Internal.DeltaDensityFitter;
import qchem.BasisSet.PointSetPullback;   // the block face leg 2 asks (lane B, READY 4)

namespace qchem::Fitting
{

DeltaDensityFitter::DeltaDensityFitter(fbs_t fbs, const BasisSet::PointCoulombIntegrals& coulomb, const BasisSet::FitQuadrature& quad)
    : itsFitBasis(std::move(fbs)), itsCoulomb(&coulomb), itsQuad(quad), itsMesh(quad.GetMesh())
{
    if (!itsFitBasis) throw std::logic_error("DeltaDensityFitter: needs a fit basis");
    if (!itsMesh) throw std::logic_error("DeltaDensityFitter: needs the fit basis's mesh");
    if (!itsFitBasis->isOrtho())
        throw std::logic_error("DeltaDensityFitter: the fit basis '"+itsFitBasis->Name()+"' is not orthogonal, so its projection "
                               "is not its fit -- a grid density fit needs a diagonal metric");
    itsW=std::make_shared<const rvec_t>(itsFitBasis->OverlapDiagonal());
    itsCharge=itsFitBasis->Charge();
    double total=0.0;
    for (size_t g=0; g<itsCharge.size(); g++) total+=itsCharge[g];
    itsC=rvec_t(itsW->size(),1.0/total);   // the initial guess: a uniform density of unit charge, as ConstrainedFF's guess is
}

double DeltaDensityFitter::HartreeEnergy(const std::function<double()>& trDV) const
{
    if (itsV.size()==0) throw std::logic_error("DeltaDensityFitter: asked for an energy before any DoFit");
    return itsCoulomb->HartreeEnergy(itsCoulomb->UsesTrDV() ? trDV() : 0.0, Projection(), itsV);
}

void DeltaDensityFitter::DoFit(const Projections& rho)
{
    itsC=rho.Vector();
    itsV=itsCoulomb->Potential(Projection());
}

hmat_t<double> DeltaDensityFitter::Repulsion(const BasisSet::Orbital_DFT_IBS<double,double>& orb) const
{
    if (itsV.size()==0) throw std::logic_error("DeltaDensityFitter: asked for a repulsion before any DoFit");
    return Integrator(orb).Adjoint(itsV);   // sum_g w_g v_g chi_i chi_j
}

bool DeltaDensityFitter::CanPullbackRepulsion(const BasisSet::Orbital_DFT_IBS<double,double>& orb) const
{
    if (!ForceGap().empty()) return false;
    const auto* face=dynamic_cast<const BasisSet::PointSetPullback<double>*>(&orb);   // abstract -> abstract, as the XC sampler's
    return face && face->CanPointSetPullback();
}

Outcome<rvec3vec_t,AD::NoGradient> DeltaDensityFitter::MakeRepulsionPullback(const BasisSet::Orbital_DFT_IBS<double,double>& orb, const hmat_t<double>& D,
                                                                              const Structure* st) const
{
    using R=Outcome<rvec3vec_t,AD::NoGradient>;
    if (const std::string gap=AtomsGap(st); !gap.empty()) return R::Fail({gap});
    const auto* face=dynamic_cast<const BasisSet::PointSetPullback<double>*>(&orb);
    if (!face || !face->CanPointSetPullback()) return R::Fail({"the orbital block "+orb.Name()+" gives no point-set pullback"});
    const qcMesh::Mesh& m=*itsMesh;
    rvec_t c(m.size());
    for (size_t g=0; g<c.size(); g++) c[g]=(*itsW)[g]*itsV[g];
    const auto p=face->MakePointSetPullback(m.Points(), c, D, st);
    if (!p) return R::Fail(p.Error());
    rvec3vec_t G=p->gradient;
    for (size_t a=0; a<m.NSites(); a++)
        for (size_t g=m.SiteBegin(a); g<m.SiteEnd(a); g++) G[a]+=p->perPoint[g];   // the point rides its owner
    ProjectOntoA1(G);   // on the representatives, c_o = W_o v_o took each orbit |o| times: g~, and P g~ is the full mesh's
    return R::Ok(std::move(G));
}

Outcome<rvec3vec_t,AD::NoGradient> DeltaDensityFitter::MakeHartreePullback(const Structure* st) const
{
    using R=Outcome<rvec3vec_t,AD::NoGradient>;
    if (const std::string gap=AtomsGap(st); !gap.empty()) return R::Fail({gap});
    const qcMesh::Mesh& m=*itsMesh;
    rvec_t f(m.size());
    for (size_t g=0; g<f.size(); g++) f[g]=itsV[g]*itsC[g];
    // Leg 1.  On the representatives: the PARENT's weight response at |o| f_o on each representative (g~), projected below.
    const qcMesh::FoldedMesh* P=itsQuad.Parent();
    const auto w= P ? P->GetMesh()->WeightResponse()->Pullback(*P->GetMesh(), P->ScatterToRepresentatives(f))
                    : m.WeightResponse()->Pullback(m, f);
    // Leg 3, scaled as qᵀKq.  A folded operator takes it on the full charges, so it is already the full mesh's (P leaves it as it is).
    const rvec3vec_t k=itsCoulomb->MakeSelfRepulsionPullback(Projection());
    rvec3vec_t G(m.NSites());
    for (size_t a=0; a<G.size(); a++) G[a]=w.gradient[a]+0.5*k[a];
    ProjectOntoA1(G);
    return R::Ok(std::move(G));
}

double DeltaDensityFitter::Integral() const
{
    double n=0.0;
    for (size_t g=0; g<itsC.size(); g++) n+=itsCharge[g]*itsC[g];
    return n;
}

std::string DeltaDensityFitter::ForceGap() const
{
    if (itsCoulomb->UsesTrDV()) return "a multipole fit's two-leg Hartree energy gives no force yet";
    if (!itsCoulomb->CanPullbackSelfRepulsion()) return "the δ fit's Coulomb operator gives no pullback (the direct sum)";
    const qcMesh::Mesh& m=*itsMesh;
    if (m.size()!=itsW->size()) return "the δ fit's functions are not its mesh's points";
    if (m.NSites()==0 || m.OwnerOffsets().size()!=m.size()) return "the δ fit's mesh records no point owners (not an atom-centred build)";
    if (const qcMesh::FoldedMesh* P=itsQuad.Parent())
    {
        const qcMesh::Mesh& pm=*P->GetMesh();
        if (!P->HasSiteAction()) return "the folded δ fit mesh's builder supplied no site action (which atom each op carries onto which)";
        if (!pm.WeightResponse()) return "the folded δ fit mesh's builder installed no weight response";
        if (pm.NSites()!=m.NSites()) return "the representatives' site blocks are not the folded mesh's";
        return "";
    }
    if (!m.WeightResponse()) return "the δ fit's mesh builder installed no weight response";
    return "";
}

std::string DeltaDensityFitter::AtomsGap(const Structure* st) const
{
    if (std::string gap=ForceGap(); !gap.empty()) return gap;
    if (itsV.size()==0) return "a δ fit asked for a force before any DoFit";
    if (!st || st->GetNumAtoms()!=itsMesh->NSites()) return "the δ fit's mesh site blocks are not the structure's atoms";
    return "";
}

const ScreenedMatrixIntegrator<double>& DeltaDensityFitter::Integrator(const BasisSet::Orbital_DFT_IBS<double,double>& orb) const
{
    std::lock_guard<std::mutex> lock(itsCacheLock.m);
    return itsInt.GetOrFill(orb.Instance(), [&]   // by the block's INSTANCE, not its symmetry (lane H 2026-11-15)
    {
        const auto& face=dynamic_cast<const BasisSet::Integrals_Overlap3C<double,double>&>(*itsFitBasis);
        return ScreenedMatrixIntegrator<double>(face.Overlap3C(orb), itsW);   // the weights are shared, not copied per block
    });
}

} //namespace
