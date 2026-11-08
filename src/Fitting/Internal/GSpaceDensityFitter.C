// File: Fitting/Internal/GSpaceDensityFitter.C  The first-order density fit on a PLANE-WAVE fit basis.
//
// V-ONE-OF-N 3b (user rulings 2026-11-08, doc/Records/VOneOfN_D2Design.md §1a, §2a row 7).  An orthonormal {G} basis is an
// ordinary fit basis whose metric happens to be diagonal, so this is the plane-wave member of the same family as the δ-mesh
// fitter (GridDensityFitter), reached through the same face (FirstOrderDensityFit) by the same term (Vee_Fit):
//   - PROJECTION.  The density projects through the OVERLAP 3-centre tensor: ρ̃(G).  On an orthonormal basis that projection IS
//     the fit (the fit metric is irrelevant, role 1 of §1a), so the coefficients are ρ̃ itself, held PRIVATELY as a ΔG_Map.
//     No kernel anywhere in this class.
//   - PHYSICS.  The repulsion is the fit basis's COULOMB 3-centre integral, ⟨ab|1/r|G⟩ = k(G)⟨ab|G⟩ (Repulsion3C's adjoint),
//     contracted with the private ρ̃ -- the physics round the Hartree term triggers by asking for it (role 2).  The energy is
//     ½ ρ̃ᵀ⟨G|1/r|G'⟩ρ̃ = ½Ω Σ k(G)|ρ̃(G)|², the basis's diagonal Coulomb metric.
//   - It replaces the hand-written V_H field the PW Hartree term used to hold (Vee_Hartree), and the kernel the density used to
//     bake into its own answer.
module;
#include <cassert>
#include <complex>
#include <functional>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <vector>
export module qchem.Fitting.Internal.GSpaceDensityFitter;
export import qchem.Fitting.CoulombOperator;   // FirstOrderDensityFit<dcmplx>
import qchem.BasisSet.Orbital_DFT_IBS;          // Orbital_DFT_IBS<U,dcmplx>, cFIT_CD_ABS / cFIT_SF_ABS
import qchem.BasisSet.PairProjection;           // GSpacePairProjectionIntegrator
import qchem.BasisSet.G_FieldEvaluator;         // G_PoissonKernel, G_RasterTransform (the PW fit basis's own faces)
import qchem.BasisSet.GMap;                     // SymmetrizeGMap (the IBZ star-average in G)
import qchem.Blaze;                             // blazem::NarrowExact
import qchem.Parallel;                          // LazyLock -- the per-block caches' guard
import qchem.Reporting;                         // report::Timed

export namespace qchem::Fitting
{

class GSpaceDensityFitter : public virtual FirstOrderDensityFit<dcmplx>
{
public:
    typedef std::shared_ptr<const BasisSet::cFIT_CD_ABS> fbs_t;
    //! The constructor DEMANDS the fit basis's four views (CLAUDE.md "Constructors pass concrete objects"): its Coulomb face
    //! \a pw (co-owned), and the overlap, Poisson-kernel and raster faces of the SAME object.  Gathering them is the
    //! factory's job (\c MakeGSpaceDensityFitter), so this class casts nothing.
    GSpaceDensityFitter(fbs_t pw, const BasisSet::cFIT_SF_ABS& sf, const BasisSet::G_PoissonKernel& pk,
                        const BasisSet::G_RasterTransform& rt)
        : itsCD(std::move(pw)), itsSF(&sf), itsPK(&pk), itsRT(&rt)
    {
        if (!itsCD) throw std::logic_error("GSpaceDensityFitter: no fit basis");
        if (!itsCD->isOrtho())
            throw std::logic_error("GSpaceDensityFitter: the fit basis '"+itsCD->Name()+"' is not orthonormal, so its projection is not its fit");
        itsVolume=itsRT->Integral(rvec_t(itsRT->RasterSize(),1.0));                 // Ω through the raster's own rule
    }

    // --- the projector the density contracts its D into (PairProjector<dcmplx>)
    const qcMesh::MatrixForward<double>& Forward(const BasisSet::Orbital_DFT_IBS<double,dcmplx>& orb) const override {return Integrator<double>(orb);}
    const qcMesh::MatrixForward<dcmplx>& Forward(const BasisSet::Orbital_DFT_IBS<dcmplx,dcmplx>& orb) const override {return Integrator<dcmplx>(orb);}
    //! A density with no D: its own Fourier components when it knows them (\c AnalyticFourier), else sampled on my raster and
    //! transformed -- either way ρ̃ on my {G}, the axis \c Forward uses.
    Projections ProjectField(const ScalarFunction<double>& f) const override
    {
        // A cross-cast between abstract faces, and the one seam there is: a density with no D reaches me as a ScalarFunction
        // through tChargeDensity::Project's default, and a second Project entry per capability is what ruling (e) declined.
        if (auto* af=dynamic_cast<const BasisSet::AnalyticFourier*>(&f)) return Projections(af->FourierComponents(*itsSF));
        return Projections(itsRT->FieldCoeffs(itsRT->ForwardFFT(itsRT->Sample(f))));
    }
    size_t NumCoefficients() const override {return itsCD->GetNumFunctions();}
    //! A zero density's ρ̃: the empty map (an absent G reads as 0).
    Projections Zero() const override {return Projections(ΔG_Map{});}
    //! The IBZ star-average of ρ̃: my functions are plane waves, so a point op permutes their G's.
    void StarAverage(Projections& p, const std::vector<Symmetry::Lattice_3D::ReciprocalOp>& ops) const override
    {
        qchem::report::Timed timed("scf: rho-tilde IBZ star-average (SymmetrizeGMap)");
        p=Projections(SymmetrizeGMap(p.GMap(), ops));
    }

    // --- the first-order fitted Hartree face
    void DoFit(const Projections& rho) override {itsRho=rho.GMap(); itsFitted=true;}
    hmat_t<double> Repulsion(const BasisSet::Orbital_DFT_IBS<double,dcmplx>& orb) const override {return RepulsionT<double>(orb);}
    hmat_t<dcmplx> Repulsion(const BasisSet::Orbital_DFT_IBS<dcmplx,dcmplx>& orb) const override {return RepulsionT<dcmplx>(orb);}
    //! \f$E_H=\tfrac12\Omega\sum_G k(G)|\tilde\rho(G)|^2\f$: first order is exact on an orthonormal basis, so \a trDV is not asked.
    double HartreeEnergy(const std::function<double()>&) const override
    {
        RequireFit("an energy");
        qchem::report::Timed timer("scf: E_H G-space pairing (no matrix)");
        double e=0.0;
        for (const auto& [dm,r] : itsRho)
            if (const double k=itsPK->CoulombKernel(dm); k>0.0) e+=k*std::norm(r);
        return 0.5*itsVolume*e;
    }
    //! \f$\int\tilde\rho=\Omega\,\tilde\rho(0)\f$.
    double Integral() const override
    {
        RequireFit("a charge");
        auto it=itsRho.find(ivec3_t(0,0,0));
        return it==itsRho.end() ? 0.0 : itsVolume*std::real(it->second);
    }

private:
    void RequireFit(const char* what) const
    { if (!itsFitted) throw std::logic_error(std::string("GSpaceDensityFitter: asked for ")+what+" before any DoFit"); }
    //! The PHYSICS round: my private ρ̃ contracted with the fit basis's Coulomb 3-centre integrals (k applied in its adjoint).
    template <class U> hmat_t<U> RepulsionT(const BasisSet::Orbital_DFT_IBS<U,dcmplx>& orb) const
    {
        RequireFit("a repulsion");
        return blazem::NarrowExact<U>(GSpacePairProjectionIntegrator<U>(orb.Repulsion3C(*itsCD)).Adjoint(itsRho));
    }
    //! The projection VIEW for one orbital block, over the basis's cached OVERLAP tensor, by spatial symmetry (as every fitter).
    template <class U> const GSpacePairProjectionIntegrator<U>& Integrator(const BasisSet::Orbital_DFT_IBS<U,dcmplx>& orb) const
    {
        auto& cache=Cache<U>();
        std::lock_guard<std::mutex> lock(itsCacheLock.m);
        const sym_t& id=orb.GetSymt();
        auto it=cache.find(id);
        if (it!=cache.end()) return it->second;
        return cache.emplace(id, GSpacePairProjectionIntegrator<U>(orb.Overlap3C(*itsSF))).first->second;
    }
    template <class U> SymMap<GSpacePairProjectionIntegrator<U>>& Cache() const
    { if constexpr (std::is_same_v<U,double>) return itsIntR; else return itsIntC; }

    fbs_t itsCD;                                    //!< the plane-wave fit basis (its Coulomb face)
    const BasisSet::cFIT_SF_ABS*       itsSF=nullptr;   //!< ...its overlap face (the projection)
    const BasisSet::G_PoissonKernel*   itsPK=nullptr;   //!< ...its diagonal Coulomb metric
    const BasisSet::G_RasterTransform* itsRT=nullptr;   //!< ...its raster
    double itsVolume=0.0;                           //!< Ω
    ΔG_Map itsRho;                                  //!< the fit = ρ̃(G), private
    bool   itsFitted=false;
    mutable SymMap<GSpacePairProjectionIntegrator<double>> itsIntR;   //!< real TRIM blocks
    mutable SymMap<GSpacePairProjectionIntegrator<dcmplx>> itsIntC;   //!< Bloch blocks
    qchem::LazyLock itsCacheLock;
};

} //namespace
