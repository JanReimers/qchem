// File: ChargeDensity/Internal/Imp/DensitySampler_OverlapProjected.C  The projected strategy over an overlap-metric fit (A1).
//
// V-XC-ONE-TERM step 3 (lane S, 2026-11-10): the XCρ::Fit cell of the molecular / atom A1 Gaussian fit, through the ONE XC term
// (tVxc_Quadrature<double>).  E = sum_g w_g e(P rho(r_g)), H = dE/dD exactly: the density's analytic projection <f_a|rho> is the
// forward, the same analytic 3-centre tensor's adjoint is H, and the mesh's Phi and w sit on both sides.  See the class comment.
module;
#include <cassert>
#include <memory>
#include <mutex>
#include <optional>
#include <stdexcept>
#include <string>
module qchem.ChargeDensity.Internal.DensitySampler;
import qchem.ChargeDensity;
import qchem.BasisSet.Orbital_DFT_IBS;
import qchem.BasisSet.IrrepBasisSet;    // Evaluatable_IBS<double>: the fit functions at points
import qchem.BasisSet.PairProjection;   // DensePairProjectionIntegrator over the basis's 3-centre overlap
import qchem.Symmetry;                  // SymMap: one integrator per orbital block, by its spatial symmetry
import qchem.Parallel;                  // LazyLock
import qchem.Projections;
import qchem.Blaze;

namespace qchem::ChargeDensity
{

void ReportNegativeRho(const DensitySamplerBase& q, const rvec_t& rho, const char* route);   // the Singles unit's

//! \brief What a density projects itself onto: \f$\langle f_a|\rho\rangle\f$ in the OVERLAP metric.  A matrix density contracts its D
//! into the per-block forward (the basis's analytic 3-centre overlap); a matrix-free one hands its field to \c ProjectField, which
//! integrates it on the XC mesh (the same \f$\Phi,w\f$ the adjoint uses).
class OverlapProjectedDensitySampler::Projector
    : public virtual BasisSet::PairProjector<double>
{
public:
    Projector(const BasisSet::rFIT_SF_ABS& fb, const rmat_t& phi, const rvec_t& w) : itsFit(fb), itsPhi(phi), itsW(w) {}
    const qcMesh::MatrixForward<double>& Forward(const BasisSet::Orbital_DFT_IBS<double,double>& orb) const override {return Integrator(orb);}
    Projections ProjectField(const ScalarFunction<double>& f) const override
    {
        const rvec3vec_t& r=Points();
        rvec_t fw=f(r);
        assert(fw.size()==itsW.size());
        for (size_t g=0; g<fw.size(); g++) fw[g]*=itsW[g];
        return Projections(rvec_t(blazem::trans(itsPhi)*fw));
    }
    size_t NumCoefficients() const override {return itsFit.GetNumFunctions();}
    //! The forward / adjoint view for one orbital block, over the basis's CACHED tensor (keyed by BasisSetID there; by the block's
    //! symmetry here, as \c ConstrainedFF does).
    const DensePairProjectionIntegrator<double>& Integrator(const BasisSet::Orbital_DFT_IBS<double,double>& orb) const
    {
        std::lock_guard<std::mutex> lock(itsLock.m);   // a parallel Fock pass's blocks share this projector
        const sym_t& id=orb.GetSymt();
        auto it=itsInt.find(id);
        if (it!=itsInt.end()) return it->second;
        return itsInt.emplace(id, DensePairProjectionIntegrator<double>(itsFit.Overlap3C(orb), itsFit.Charge())).first->second;
    }
    void SetPoints(const rvec3vec_t* r) {itsPoints=r;}
private:
    const rvec3vec_t& Points() const {assert(itsPoints); return *itsPoints;}
    const BasisSet::rFIT_SF_ABS& itsFit;
    const rmat_t&                itsPhi;
    const rvec_t&                itsW;
    const rvec3vec_t*            itsPoints=nullptr;
    mutable SymMap<DensePairProjectionIntegrator<double>> itsInt;
    qchem::LazyLock itsLock;
};

OverlapProjectedDensitySampler::OverlapProjectedDensitySampler(fbs_t fb, BasisSet::FitQuadrature quad)
    : itsFitBasis(std::move(fb))
    , itsMetric(dynamic_cast<const BasisSet::FIT_SF_NonOrtho*>(itsFitBasis.get()))   // abstract -> abstract
    , itsQuad(std::move(quad))
{
    assert(itsFitBasis);
    if (!itsMetric)
        throw std::logic_error("OverlapProjectedDensitySampler: the fit basis has no overlap metric (FIT_SF_NonOrtho) to project onto");
    const auto* eval=dynamic_cast<const BasisSet::Evaluatable_IBS<double>*>(itsFitBasis.get());   // abstract -> abstract
    if (!eval)
        throw std::logic_error("OverlapProjectedDensitySampler: the fit basis cannot evaluate its functions at points (Evaluatable_IBS), "
                               "so P rho has no values on the XC mesh");
    if (!itsQuad.GetMesh())
        throw std::logic_error("OverlapProjectedDensitySampler: no XC mesh was handed in -- the projected strategy integrates on it");
    const qcMesh::Mesh& m=*itsQuad.GetMesh();
    itsPhi=(*eval)(m.Points());   // points x functions
    if (itsPhi.columns()!=itsFitBasis->GetNumFunctions())
        throw std::logic_error("OverlapProjectedDensitySampler: the fit basis evaluated a different number of functions than it holds");
    itsProjector=std::make_unique<Projector>(*itsFitBasis, itsPhi, m.Weights());
    itsProjector->SetPoints(&m.Points());
}
OverlapProjectedDensitySampler::~OverlapProjectedDensitySampler() = default;

double OverlapProjectedDensitySampler::Integrate(const rvec_t& f) const
{
    const rvec_t& w=itsQuad.GetMesh()->Weights();
    assert(f.size()==w.size());
    return blazem::dot(w, f);
}
size_t OverlapProjectedDensitySampler::NumPoints() const {return itsQuad.GetMesh()->size();}

rvec_t OverlapProjectedDensitySampler::SiteIntegrals(const rvec_t& f) const
{
    const qcMesh::Mesh& m=*itsQuad.GetMesh();
    rvec_t ret(m.NSites(), 0.0);
    const rvec_t& w=m.Weights();
    for (size_t a=0; a<m.NSites(); a++)
        for (size_t g=m.SiteBegin(a); g<(a+1<m.NSites() ? m.SiteBegin(a+1) : m.size()); g++) ret[a]+=w[g]*f[g];
    return ret;
}

rvec_t OverlapProjectedDensitySampler::SampleOne(const rChargeDensity* cd) const
{
    assert(cd);
    const Projections p=cd->Project(*itsProjector);
    const rvec_t c=itsMetric->InvOverlap()*p.Vector();
    return itsPhi*c;
}

void OverlapProjectedDensitySampler::Refresh(const rChargeDensity* cd) const
{
    assert(itsPolVersion==size_t(-1) && "OverlapProjectedDensitySampler: this engine already served RhoPol (no cross-invalidation)");
    if (cd->Version()==itsRhoVersion) return;
    itsRhoVersion=cd->Version();
    itsRho=SampleOne(cd);
    ReportNegativeRho(*this, itsRho, "overlap-projected");
}

void OverlapProjectedDensitySampler::RefreshPol(const rChargeDensity* cd) const
{
    assert(itsRhoVersion==size_t(-1) && "OverlapProjectedDensitySampler: this engine already served Rho (no cross-invalidation)");
    if (cd->Version()==itsPolVersion) return;
    itsPolVersion=cd->Version();
    const rChargeDensity* up=ChannelOf(cd, Spin::Up);
    const rChargeDensity* dn=ChannelOf(cd, Spin::Down);
    if (!up || !dn)
    {
        // The spin-agnostic seed: rho_up = rho_dn = rho/2 (P is linear, so P(rho/2) = (P rho)/2).
        itsRhoUp=SampleOne(cd);
        itsRhoUp*=0.5;
        itsRhoDn=itsRhoUp;
        ReportNegativeRho(*this, itsRhoUp, "overlap-projected half-density seed");
        return;
    }
    itsRhoUp=SampleOne(up);
    itsRhoDn=SampleOne(dn);
    ReportNegativeRho(*this, itsRhoUp, "overlap-projected(up)");
    ReportNegativeRho(*this, itsRhoDn, "overlap-projected(dn)");
}

const rvec_t& OverlapProjectedDensitySampler::Rho(const rChargeDensity* cd) const {Refresh(cd); return itsRho;}
const rvec_t& OverlapProjectedDensitySampler::RhoPol(const rChargeDensity* cd, const Spin& s) const
{
    assert(s!=Spin::None && "OverlapProjectedDensitySampler::RhoPol: ask for a channel, not the total");
    RefreshPol(cd);
    return s==Spin::Up ? itsRhoUp : itsRhoDn;
}

// H = sum_a (S^-1 Phi^T (w v))_a <ij|a>: the exact adjoint of SampleOne's P rho on the same points and weights.
rsmat_t OverlapProjectedDensitySampler::Matrix(const robs_t* bs, const rvec_t& v) const
{
    const auto& orb=dynamic_cast<const BasisSet::Orbital_DFT_IBS<double,double>&>(*bs);   // genuine "is it?" cross-cast (throws)
    const rvec_t& w=itsQuad.GetMesh()->Weights();
    assert(v.size()==w.size());
    rvec_t wv(v.size());
    for (size_t g=0; g<v.size(); g++) wv[g]=w[g]*v[g];
    const rvec_t cv=itsMetric->InvOverlap()*rvec_t(blazem::trans(itsPhi)*wv);
    const rsmat_t H=itsProjector->Integrator(orb).Adjoint(cv);
    assert(!blazem::isnan(H));
    return H;
}

} //namespace
