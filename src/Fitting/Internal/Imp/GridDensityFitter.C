// File: Fitting/Internal/Imp/GridDensityFitter.C  The bodies of qchem.Fitting.Internal.GridDensityFitter: a density fit on a mesh, its
// Hartree energy and its fixed-D force legs (the interface carries the design and the rulings).
module;
#include <functional>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
module qchem.Fitting.Internal.GridDensityFitter;
import qchem.BasisSet.PointSetPullback;   // the block face leg 2 asks (lane B, READY 4)

namespace qchem::Fitting
{

GridDensityFitter::GridDensityFitter(fbs_t fbs, const BasisSet::PointCoulombIntegrals& coulomb, std::shared_ptr<const qcMesh::Mesh> mesh)
    : itsFitBasis(std::move(fbs)), itsCoulomb(&coulomb), itsMesh(std::move(mesh))
{
    if (!itsFitBasis) throw std::logic_error("GridDensityFitter: needs a fit basis");
    if (!itsMesh) throw std::logic_error("GridDensityFitter: needs the fit basis's mesh");
    if (!itsFitBasis->isOrtho())
        throw std::logic_error("GridDensityFitter: the fit basis '"+itsFitBasis->Name()+"' is not orthogonal, so its projection "
                               "is not its fit -- a grid density fit needs a diagonal metric");
    itsW=std::make_shared<const rvec_t>(itsFitBasis->OverlapDiagonal());
    itsCharge=itsFitBasis->Charge();
    double total=0.0;
    for (size_t g=0; g<itsCharge.size(); g++) total+=itsCharge[g];
    itsC=rvec_t(itsW->size(),1.0/total);   // the initial guess: a uniform density of unit charge, as ConstrainedFF's guess is
}

double GridDensityFitter::HartreeEnergy(const std::function<double()>& trDV) const
{
    if (itsV.size()==0) throw std::logic_error("GridDensityFitter: asked for an energy before any DoFit");
    return itsCoulomb->HartreeEnergy(itsCoulomb->UsesTrDV() ? trDV() : 0.0, Projection(), itsV);
}

void GridDensityFitter::DoFit(const Projections& rho)
{
    itsC=rho.Vector();
    itsV=itsCoulomb->Potential(Projection());
}

hmat_t<double> GridDensityFitter::Repulsion(const BasisSet::Orbital_DFT_IBS<double,double>& orb) const
{
    if (itsV.size()==0) throw std::logic_error("GridDensityFitter: asked for a repulsion before any DoFit");
    return Integrator(orb).Adjoint(itsV);   // sum_g w_g v_g chi_i chi_j
}

bool GridDensityFitter::CanPullbackRepulsion(const BasisSet::Orbital_DFT_IBS<double,double>& orb) const
{
    if (!ForceGap().empty()) return false;
    const auto* face=dynamic_cast<const BasisSet::PointSetPullback<double>*>(&orb);   // abstract -> abstract, as the XC sampler's
    return face && face->CanPointSetPullback();
}

Outcome<rvec3vec_t,AD::NoGradient> GridDensityFitter::MakeRepulsionPullback(const BasisSet::Orbital_DFT_IBS<double,double>& orb, const hmat_t<double>& D,
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
    return R::Ok(std::move(G));
}

Outcome<rvec3vec_t,AD::NoGradient> GridDensityFitter::MakeHartreePullback(const Structure* st) const
{
    using R=Outcome<rvec3vec_t,AD::NoGradient>;
    if (const std::string gap=AtomsGap(st); !gap.empty()) return R::Fail({gap});
    const qcMesh::Mesh& m=*itsMesh;
    rvec_t f(m.size());
    for (size_t g=0; g<f.size(); g++) f[g]=itsV[g]*itsC[g];
    const auto w=m.WeightResponse()->Pullback(m, f);                       // leg 1
    const rvec3vec_t k=itsCoulomb->MakeSelfRepulsionPullback(Projection());  // leg 3, scaled as qᵀKq
    rvec3vec_t G(m.NSites());
    for (size_t a=0; a<G.size(); a++) G[a]=w.gradient[a]+0.5*k[a];
    return R::Ok(std::move(G));
}

double GridDensityFitter::Integral() const
{
    double n=0.0;
    for (size_t g=0; g<itsC.size(); g++) n+=itsCharge[g]*itsC[g];
    return n;
}

std::string GridDensityFitter::ForceGap() const
{
    if (itsCoulomb->UsesTrDV()) return "a multipole fit's two-leg Hartree energy gives no force yet";
    if (!itsCoulomb->CanPullbackSelfRepulsion()) return "the δ fit's Coulomb operator gives no pullback (the direct sum)";
    const qcMesh::Mesh& m=*itsMesh;
    if (m.size()!=itsW->size()) return "the δ fit's functions are not its mesh's points";
    if (m.NSites()==0 || m.OwnerOffsets().size()!=m.size()) return "the δ fit's mesh records no point owners (not an atom-centred build)";
    if (!m.WeightResponse()) return "the δ fit's mesh builder installed no weight response";
    return "";
}

std::string GridDensityFitter::AtomsGap(const Structure* st) const
{
    if (std::string gap=ForceGap(); !gap.empty()) return gap;
    if (itsV.size()==0) return "a δ fit asked for a force before any DoFit";
    if (!st || st->GetNumAtoms()!=itsMesh->NSites()) return "the δ fit's mesh site blocks are not the structure's atoms";
    return "";
}

const ScreenedMatrixIntegrator<double>& GridDensityFitter::Integrator(const BasisSet::Orbital_DFT_IBS<double,double>& orb) const
{
    std::lock_guard<std::mutex> lock(itsCacheLock.m);
    return itsInt.GetOrFill(orb.Instance(), [&]   // by the block's INSTANCE, not its symmetry (lane H 2026-11-15)
    {
        const auto& face=dynamic_cast<const BasisSet::Integrals_Overlap3C<double,double>&>(*itsFitBasis);
        return ScreenedMatrixIntegrator<double>(face.Overlap3C(orb), itsW);   // the weights are shared, not copied per block
    });
}

} //namespace
