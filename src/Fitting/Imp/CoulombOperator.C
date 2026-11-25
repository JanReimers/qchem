// File: Fitting/Imp/CoulombOperator.C  The factories of the grid density route.
module;
#include <memory>
#include <vector>
#include <stdexcept>
#include <string>
module qchem.Fitting.CoulombOperator;
import qchem.Fitting.Internal.DeltaDensityFitter;
import qchem.Fitting.Internal.DeltaProjectionDensityFitter;
import qchem.Fitting.Internal.GSpaceDensityFitter;
import qchem.Fitting.Internal.FunctionFitterImp;   // IntegralConstrainedFF (the Gaussian Dunlap fit)
import qchem.BasisSet.G_FieldEvaluator;
import qchem.BasisSet.GSpaceDensityProjector;   // the ρ̃ projector factory (the fitter's projection half)
import qchem.BasisSet.CoulombDeltaFit_IBS;   // PointCoulombIntegrals   // the plane-wave fit basis faces the factory gathers

namespace qchem::Fitting
{

// The whole-density defaults (the k-fold's per-block fallback): each block through its own derivative, summed.
Outcome<rvec3vec_t,AD::NoGradient> FirstOrderDensityFit<dcmplx>::MakeRepulsionPullbackBlocks(const std::vector<BasisSet::FoldBlockRef>& blocks,
                                                                                            const std::vector<BasisSet::FoldBlockD>& Ds, const Structure* st) const
{
    using R=Outcome<rvec3vec_t,AD::NoGradient>;
    if (blocks.size()!=Ds.size()) throw std::logic_error("MakeRepulsionPullbackBlocks: "+std::to_string(Ds.size())+" density matrices for "+std::to_string(blocks.size())+" blocks");
    rvec3vec_t F;
    for (size_t b=0;b<blocks.size();b++)
    {
        const R r= (blocks[b].real && Ds[b].real)   ? MakeRepulsionPullback(*blocks[b].real,  *Ds[b].real,  st)
                 : (blocks[b].bloch && Ds[b].bloch) ? MakeRepulsionPullback(*blocks[b].bloch, *Ds[b].bloch, st)
                 : throw std::logic_error("MakeRepulsionPullbackBlocks: block "+std::to_string(b)+"'s density matrix is not on the block's own face");
        if (!r) return r;
        if (F.size()==0) { F=*r; continue; }
        if (r->size()!=F.size()) throw std::logic_error("MakeRepulsionPullbackBlocks: two blocks answered for different atom counts");
        for (size_t A=0;A<F.size();A++) F[A]+=(*r)[A];
    }
    return R::Ok(std::move(F));
}
Outcome<rmat3d_t,AD::NoGradient> FirstOrderDensityFit<dcmplx>::MakeRepulsionVirialBlocks(const std::vector<BasisSet::FoldBlockRef>& blocks,
                                                                                        const std::vector<BasisSet::FoldBlockD>& Ds, const Structure* st) const
{
    using R=Outcome<rmat3d_t,AD::NoGradient>;
    if (blocks.size()!=Ds.size()) throw std::logic_error("MakeRepulsionVirialBlocks: "+std::to_string(Ds.size())+" density matrices for "+std::to_string(blocks.size())+" blocks");
    rmat3d_t V(0,0,0, 0,0,0, 0,0,0);
    for (size_t b=0;b<blocks.size();b++)
    {
        const R r= (blocks[b].real && Ds[b].real)   ? MakeRepulsionVirial(*blocks[b].real,  *Ds[b].real,  st)
                 : (blocks[b].bloch && Ds[b].bloch) ? MakeRepulsionVirial(*blocks[b].bloch, *Ds[b].bloch, st)
                 : throw std::logic_error("MakeRepulsionVirialBlocks: block "+std::to_string(b)+"'s density matrix is not on the block's own face");
        if (!r) return r;
        V+=*r;
    }
    return R::Ok(V);
}

std::unique_ptr<FirstOrderDensityFit<double>>
MakeDeltaDensityFitter(std::shared_ptr<const BasisSet::rFIT_SF_ABS> deltaBasis, std::shared_ptr<const qcMesh::Mesh> mesh)
{
    // The factory gathers the basis's Coulomb-integral face (the concrete association lives here, not in the class).
    if (!deltaBasis) throw std::logic_error("MakeDeltaDensityFitter: no fit basis");
    auto* coulomb=dynamic_cast<const BasisSet::PointCoulombIntegrals*>(deltaBasis.get());
    if (!coulomb)
        throw std::logic_error("MakeDeltaDensityFitter: the fit basis '"+deltaBasis->Name()+"' carries no Coulomb integrals "
                               "(wrap a δ basis in BasisSet::CoulombDeltaFit_IBS with its operator)");
    if (!mesh) throw std::logic_error("MakeDeltaDensityFitter: no mesh for the fit basis '"+deltaBasis->Name()+"'");
    return std::make_unique<DeltaDensityFitter>(std::move(deltaBasis), *coulomb, std::move(mesh));
}

std::unique_ptr<FirstOrderDensityFit<double>>
MakeDeltaProjectionDensityFitter(std::shared_ptr<const BasisSet::rFIT_SF_ABS> deltaBasis, std::shared_ptr<const qcMesh::Mesh> mesh)
{
    if (!deltaBasis) throw std::logic_error("MakeDeltaProjectionDensityFitter: no fit basis");
    if (!mesh) throw std::logic_error("MakeDeltaProjectionDensityFitter: no mesh for the fit basis '"+deltaBasis->Name()+"'");
    return std::make_unique<DeltaProjectionDensityFitter>(std::move(deltaBasis), std::move(mesh));
}

std::unique_ptr<DunlapDensityFit<double>> MakeDunlapDensityFitter(std::shared_ptr<const BasisSet::rFIT_CD_ABS> fbs, DensityFitConstraint constraint)
{
    if (!fbs) throw std::logic_error("MakeDunlapDensityFitter: no fit basis");
    auto nonOrtho=std::dynamic_pointer_cast<const BasisSet::FIT_CD_NonOrtho>(fbs);   // its Coulomb metric-solve face
    if (fbs->isOrtho() || !nonOrtho)
        throw std::logic_error("MakeDunlapDensityFitter: the fit basis '"+fbs->Name()+"' has no Coulomb metric to solve in");
    return std::make_unique<IntegralConstrainedFF<double>>(nonOrtho, constraint==DensityFitConstraint::Charge);
}

std::unique_ptr<FirstOrderDensityFit<dcmplx>> MakeGSpaceDensityFitter(std::shared_ptr<const BasisSet::cFIT_CD_ABS> pw)
{
    // The factory gathers the plane-wave fit basis's faces (the concrete association lives here, not in the class).
    if (!pw) throw std::logic_error("MakeGSpaceDensityFitter: no fit basis");
    auto* sf=dynamic_cast<const BasisSet::cFIT_SF_ABS*>(pw.get());
    auto* pk=dynamic_cast<const BasisSet::G_PoissonKernel*>(pw.get());
    auto* rt=dynamic_cast<const BasisSet::G_RasterTransform*>(pw.get());
    if (!sf || !pk || !rt)
        throw std::logic_error("MakeGSpaceDensityFitter: the fit basis '"+pw->Name()+"' lacks a G-space face (overlap, Poisson kernel, raster)");
    // The projection half: ρ̃ on the same object's overlap face, co-owned through pw (an aliasing shared_ptr).
    auto projector=BasisSet::MakeGSpaceDensityProjector(std::shared_ptr<const BasisSet::cFIT_SF_ABS>(pw, sf),
                                                       "scf: rho-tilde IBZ star-average (SymmetrizeGMap)");
    return std::make_unique<GSpaceDensityFitter>(std::move(pw), *pk, *rt, std::move(projector));
}

} //namespace
