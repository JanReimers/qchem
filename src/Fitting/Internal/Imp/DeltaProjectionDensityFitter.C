// File: Fitting/Internal/Imp/DeltaProjectionDensityFitter.C  The bodies of qchem.Fitting.Internal.DeltaProjectionDensityFitter (the interface
// carries the derivation and the ruling).
module;
#include <functional>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
module qchem.Fitting.Internal.DeltaProjectionDensityFitter;

namespace qchem::Fitting
{

//! One orbital block: its overlap projection (the δ basis's cached \f$\langle ij|\delta_g\rangle\f$), its Coulomb projection (the block's
//! \f$(ij|\delta_g)\f$, matrix-free, owned here), and the joint forward the density contracts its D into.
struct DeltaProjectionDensityFitter::Block
    : public virtual qcMesh::MatrixForward<double>
{
    Block(const PairProjection<double>& overlap, PairProjection<double> coulomb, std::shared_ptr<const rvec_t> w)
        : coul(std::make_unique<const PairProjection<double>>(std::move(coulomb)))
        , ovl(overlap, w), cou(*coul, w) {}

    using qcMesh::MatrixForward<double>::Forward;   // un-hide the factored overload (qchem.Mesh.Integrator's warning)
    //! \f$[\rho(r_g);\ \tilde v_g]\f$
    Projections Forward(const hmat_t<double>& D) const override
    {
        const rvec_t rho=ovl.Forward(D).Vector(), v=cou.Forward(D).Vector();
        const size_t n=rho.size();
        rvec_t out(2*n);
        for (size_t g=0;g<n;g++) { out[g]=rho[g]; out[n+g]=v[g]; }
        return Projections(std::move(out));
    }
    size_t NumCoefficients() const override {return 2*ovl.NumCoefficients();}   //!< [ρ ; ṽ]

    std::unique_ptr<const PairProjection<double>> coul;   //!< owned: the Coulomb integrator borrows it
    ScreenedMatrixIntegrator<double> ovl, cou;
};

DeltaProjectionDensityFitter::DeltaProjectionDensityFitter(fbs_t fbs, std::shared_ptr<const qcMesh::Mesh> mesh)
    : itsFitBasis(std::move(fbs)), itsMesh(std::move(mesh))
{
    if (!itsFitBasis) throw std::logic_error("DeltaProjectionDensityFitter: needs a fit basis");
    if (!itsMesh) throw std::logic_error("DeltaProjectionDensityFitter: needs the fit basis's mesh");
    if (!itsFitBasis->isOrtho())
        throw std::logic_error("DeltaProjectionDensityFitter: the fit basis '"+itsFitBasis->Name()+"' is not orthogonal, so its projection "
                               "is not its fit -- a δ-projection fit needs a diagonal metric");
    itsW=std::make_shared<const rvec_t>(itsFitBasis->OverlapDiagonal());
    if (itsMesh->size()!=itsW->size())
        throw std::logic_error("DeltaProjectionDensityFitter: the mesh has "+std::to_string(itsMesh->size())+" points for the fit basis's "
                               +std::to_string(itsW->size())+" functions -- the Coulomb integrals must be taken at the fit's own points");
    itsCharge=itsFitBasis->Charge();
}

DeltaProjectionDensityFitter::~DeltaProjectionDensityFitter() = default;

double DeltaProjectionDensityFitter::HartreeEnergy(const std::function<double()>&) const
{
    if (itsV.size()==0) throw std::logic_error("DeltaProjectionDensityFitter: asked for an energy before any DoFit");
    double e=0.0;
    for (size_t g=0;g<itsQ.size();g++) e+=itsQ[g]*itsV[g];
    return 0.5*e;
}

void DeltaProjectionDensityFitter::DoFit(const Projections& proj)
{
    const rvec_t& p=proj.Vector();
    const size_t n=itsW->size();
    if (p.size()!=2*n)
        throw std::logic_error("DeltaProjectionDensityFitter::DoFit: "+std::to_string(p.size())+" projections for "+std::to_string(n)
                               +" points -- expected [rho ; v] (two per point, from this fitter's own Forward)");
    itsRho=rvec_t(n); itsQ=rvec_t(n); itsV=rvec_t(n);
    for (size_t g=0;g<n;g++) { itsRho[g]=p[g]; itsQ[g]=(*itsW)[g]*p[g]; itsV[g]=p[n+g]; }
}

hmat_t<double> DeltaProjectionDensityFitter::Repulsion(const BasisSet::Orbital_DFT_IBS<double,double>& orb) const
{
    if (itsV.size()==0) throw std::logic_error("DeltaProjectionDensityFitter: asked for a repulsion before any DoFit");
    const Block& b=BlockFor(orb);
    // ½ Σ_g (w_g ṽ_g) ⟨ij|δ_g⟩: the δ basis's overlap adjoint applies the weights itself (H = Σ_g w_g c_g χ_iχ_j(r_g), DeltaFit_IBS), so it
    // takes c = ½ṽ; the Coulomb adjoint is the bare Σ_g c_g (ij|δ_g), so it takes c = ½q.
    rvec_t c(itsV.size()), q2(itsQ.size());
    for (size_t g=0;g<c.size();g++) { c[g]=0.5*itsV[g]; q2[g]=0.5*itsQ[g]; }
    return b.ovl.Adjoint(c) + b.cou.Adjoint(q2);   // + ½ Σ_g q_g (ij|δ_g)
}

Outcome<rvec3vec_t,AD::NoGradient> DeltaProjectionDensityFitter::MakeRepulsionPullback(const BasisSet::Orbital_DFT_IBS<double,double>&,
                                                                                        const hmat_t<double>&, const Structure*) const
{
    return Outcome<rvec3vec_t,AD::NoGradient>::Fail({"the δ-projection Hartree's fixed-D force is not built yet (lane-s-design-ij-delta.md §4)"});
}

double DeltaProjectionDensityFitter::Integral() const
{
    double n=0.0;
    for (size_t g=0;g<itsRho.size();g++) n+=itsCharge[g]*itsRho[g];   // Σ ⟨δ_g|1⟩ ρ(r_g)
    return n;
}

const qcMesh::MatrixForward<double>& DeltaProjectionDensityFitter::Forward(const BasisSet::Orbital_DFT_IBS<double,double>& orb) const
{
    return BlockFor(orb);
}

Projections DeltaProjectionDensityFitter::ProjectField(const ScalarFunction<double>&) const
{
    throw std::logic_error("DeltaProjectionDensityFitter: a density with no density matrix (a field) has no analytic (rho|delta_g); "
                           "the δ-projection Hartree needs D");
}

const DeltaProjectionDensityFitter::Block& DeltaProjectionDensityFitter::BlockFor(const BasisSet::Orbital_DFT_IBS<double,double>& orb) const
{
    std::lock_guard<std::mutex> lock(itsCacheLock.m);
    return *itsBlocks.GetOrFill(orb.Instance(), [&]   // by the block's INSTANCE, as DeltaDensityFitter
    {
        const auto& ovlFace=dynamic_cast<const BasisSet::Integrals_Overlap3C<double,double>&>(*itsFitBasis);
        const auto* pc=dynamic_cast<const BasisSet::PointCoulomb<double>*>(&orb);   // abstract -> abstract, as the XC sampler's
        if (!pc || !pc->CanPointCoulomb())
            throw std::logic_error("DeltaProjectionDensityFitter: the orbital block "+orb.Name()+" gives no point-charge Coulomb integrals "
                                   "(BasisSet::PointCoulomb), so the δ-projection Hartree cannot run on it");
        return std::make_unique<Block>(ovlFace.Overlap3C(orb), pc->MakePointCoulomb(itsMesh->Points()), itsW);
    });
}

} //namespace
