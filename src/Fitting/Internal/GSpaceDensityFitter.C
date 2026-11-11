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
import qchem.Blaze;                             // blazem::NarrowExact
import qchem.Parallel;                          // LazyLock -- the per-block caches' guard
import qchem.Reporting;                         // report::Timed
import qchem.Structure;                         // the pullbacks' atoms

export namespace qchem::Fitting
{

class GSpaceDensityFitter : public virtual FirstOrderDensityFit<dcmplx>
{
public:
    typedef std::shared_ptr<const BasisSet::cFIT_CD_ABS> fbs_t;
    typedef std::unique_ptr<const BasisSet::PairProjector<dcmplx>> projector_t;
    //! The constructor DEMANDS what it uses (CLAUDE.md "Constructors pass concrete objects"): the fit basis's Coulomb face \a pw
    //! (co-owned), its Poisson-kernel and raster faces of the SAME object, and the ρ̃ projector \a projector over that basis's
    //! overlap face (the projection half; \c BasisSet::MakeGSpaceDensityProjector).  Gathering them is the factory's job
    //! (\c MakeGSpaceDensityFitter), so this class casts nothing.
    GSpaceDensityFitter(fbs_t pw, const BasisSet::G_PoissonKernel& pk, const BasisSet::G_RasterTransform& rt, projector_t projector)
        : itsCD(std::move(pw)), itsPK(&pk), itsRT(&rt), itsProjector(std::move(projector))
    {
        if (!itsProjector) throw std::logic_error("GSpaceDensityFitter: no projector");
        if (!itsCD) throw std::logic_error("GSpaceDensityFitter: no fit basis");
        if (!itsCD->isOrtho())
            throw std::logic_error("GSpaceDensityFitter: the fit basis '"+itsCD->Name()+"' is not orthonormal, so its projection is not its fit");
        itsVolume=itsRT->Integral(rvec_t(itsRT->RasterSize(),1.0));                 // Ω through the raster's own rule
    }

    // --- the projector the density contracts its D into (PairProjector<dcmplx>): ρ̃ on my {G}, by the shared projector
    const qcMesh::MatrixForward<double>& Forward(const BasisSet::Orbital_DFT_IBS<double,dcmplx>& orb) const override {return itsProjector->Forward(orb);}
    const qcMesh::MatrixForward<dcmplx>& Forward(const BasisSet::Orbital_DFT_IBS<dcmplx,dcmplx>& orb) const override {return itsProjector->Forward(orb);}
    Projections ProjectField(const ScalarFunction<double>& f) const override {return itsProjector->ProjectField(f);}
    size_t NumCoefficients() const override {return itsProjector->NumCoefficients();}
    Projections Zero() const override {return itsProjector->Zero();}
    void StarAverage(Projections& p, const Symmetry::StarGroup& group) const override
        {itsProjector->StarAverage(p,group);}

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

    //! The pullback of the physics round at my FIXED ρ̃: the fit functions are plane waves (no centre), so only the block's moving
    //! centres contribute, through the Coulomb tensor's own derivative realization (GPW: the jet walk of the gather; a plane-wave
    //! block: exactly zero).  A real TRIM block's D is promoted: the tensor is complex for both block scalars.
    bool CanPullbackRepulsion(const BasisSet::Orbital_DFT_IBS<double,dcmplx>& orb) const override {return RealisesAdjointPullback(orb.Repulsion3C(*itsCD));}
    bool CanPullbackRepulsion(const BasisSet::Orbital_DFT_IBS<dcmplx,dcmplx>& orb) const override {return RealisesAdjointPullback(orb.Repulsion3C(*itsCD));}
    Outcome<rvec3vec_t,AD::NoGradient> MakeRepulsionPullback(const BasisSet::Orbital_DFT_IBS<double,dcmplx>& orb, const hmat_t<double>& D,
                                                             const Structure* st) const override
    {
        chmat_t Dc(D.rows());
        for (size_t i=0;i<D.rows();i++) for (size_t j=i;j<D.rows();j++) Dc(i,j)=dcmplx(D(i,j),0.0);
        return RepulsionPullbackT(orb, Dc, st);
    }
    Outcome<rvec3vec_t,AD::NoGradient> MakeRepulsionPullback(const BasisSet::Orbital_DFT_IBS<dcmplx,dcmplx>& orb, const hmat_t<dcmplx>& D,
                                                             const Structure* st) const override
    {return RepulsionPullbackT(orb, D, st);}

    //! THE FIT SIDE of the virial (CellHalf_Design §5): \f$E_H=\tfrac1{2\Omega}\sum_G k(G)|N_G|^2\f$ with the counts \f$N_G=\Omega\tilde\rho(G)\f$
    //! held, so \f$\partial E_H/\partial\varepsilon_{ab}=-\delta_{ab}E_H+\tfrac12\Omega\sum_G\partial_{ab}k(G)\,|\tilde\rho(G)|^2\f$.
    bool CanHartreeVirial() const override {return true;}
    Outcome<rmat3d_t,AD::NoGradient> MakeHartreeVirial() const override
    {
        RequireFit("a Hartree virial");
        const double E=HartreeEnergy({});
        rmat3d_t V(0,0,0, 0,0,0, 0,0,0);
        for (const auto& [dm,r] : itsRho)
        {
            const rmat3d_t dk=itsPK->CoulombKernelStrain(dm);
            const double w=0.5*itsVolume*std::norm(r);
            for (size_t a=1;a<=3;a++) for (size_t b=1;b<=3;b++) V(a,b)+=w*dk(a,b);
        }
        for (size_t a=1;a<=3;a++) V(a,a)-=E;
        return Outcome<rmat3d_t,AD::NoGradient>::Ok(V);
    }
    //! THE ORBITAL SIDE (CellHalf_Design §5): the Coulomb tensor's own realization at my fixed ρ̃ (a plane-wave block: exactly zero).
    bool CanVirialRepulsion(const BasisSet::Orbital_DFT_IBS<double,dcmplx>& orb) const override {return RealisesAdjointVirial(orb.Repulsion3C(*itsCD));}
    bool CanVirialRepulsion(const BasisSet::Orbital_DFT_IBS<dcmplx,dcmplx>& orb) const override {return RealisesAdjointVirial(orb.Repulsion3C(*itsCD));}
    Outcome<rmat3d_t,AD::NoGradient> MakeRepulsionVirial(const BasisSet::Orbital_DFT_IBS<double,dcmplx>& orb, const hmat_t<double>& D,
                                                         const Structure* st) const override
    {
        chmat_t Dc(D.rows());
        for (size_t i=0;i<D.rows();i++) for (size_t j=i;j<D.rows();j++) Dc(i,j)=dcmplx(D(i,j),0.0);
        return RepulsionVirialT(orb, Dc, st);
    }
    Outcome<rmat3d_t,AD::NoGradient> MakeRepulsionVirial(const BasisSet::Orbital_DFT_IBS<dcmplx,dcmplx>& orb, const hmat_t<dcmplx>& D,
                                                         const Structure* st) const override
    {return RepulsionVirialT(orb, D, st);}

private:
    template <class U> Outcome<rmat3d_t,AD::NoGradient> RepulsionVirialT(const BasisSet::Orbital_DFT_IBS<U,dcmplx>& orb, const chmat_t& D,
                                                                         const Structure* st) const
    {
        RequireFit("a repulsion virial");
        const PairProjection<dcmplx>& g=orb.Repulsion3C(*itsCD);
        if (!RealisesAdjointVirial(g))
            return Outcome<rmat3d_t,AD::NoGradient>::Fail({"the block "+orb.Name()+"'s Coulomb tensor against the plane-wave fit gives no virial"});
        const ΔG_Map& rho=itsRho;
        return Outcome<rmat3d_t,AD::NoGradient>::Ok(ContractAdjointVirial(g, [&rho](const ivec3_t& dm)->dcmplx
            { auto it=rho.find(dm); return it==rho.end() ? dcmplx(0.0) : it->second; }, D, st));
    }
    template <class U> Outcome<rvec3vec_t,AD::NoGradient> RepulsionPullbackT(const BasisSet::Orbital_DFT_IBS<U,dcmplx>& orb, const chmat_t& D,
                                                                             const Structure* st) const
    {
        RequireFit("a repulsion pullback");
        const PairProjection<dcmplx>& g=orb.Repulsion3C(*itsCD);
        if (!RealisesAdjointPullback(g))
            return Outcome<rvec3vec_t,AD::NoGradient>::Fail({"the block "+orb.Name()+"'s Coulomb tensor against the plane-wave fit gives no pullback"});
        const ΔG_Map& rho=itsRho;
        return Outcome<rvec3vec_t,AD::NoGradient>::Ok(ContractAdjointPullback(g, [&rho](const ivec3_t& dm)->dcmplx
            { auto it=rho.find(dm); return it==rho.end() ? dcmplx(0.0) : it->second; }, D, st));
    }
    void RequireFit(const char* what) const
    { if (!itsFitted) throw std::logic_error(std::string("GSpaceDensityFitter: asked for ")+what+" before any DoFit"); }
    //! The PHYSICS round: my private ρ̃ contracted with the fit basis's Coulomb 3-centre integrals (k applied in its adjoint).
    template <class U> hmat_t<U> RepulsionT(const BasisSet::Orbital_DFT_IBS<U,dcmplx>& orb) const
    {
        RequireFit("a repulsion");
        // Adjoint-only: the Coulomb tensor has no forward (V-ONE-OF-N 3c), so no integrator pair is built over it.
        const ΔG_Map& rho=itsRho;
        return blazem::NarrowExact<U>(ContractAdjoint(orb.Repulsion3C(*itsCD), [&rho](const ivec3_t& dm)->dcmplx
            { auto it=rho.find(dm); return it==rho.end() ? dcmplx(0.0) : it->second; }));
    }
    fbs_t itsCD;                                    //!< the plane-wave fit basis (its Coulomb face)
    const BasisSet::G_PoissonKernel*   itsPK=nullptr;   //!< ...its diagonal Coulomb metric
    const BasisSet::G_RasterTransform* itsRT=nullptr;   //!< ...its raster
    double itsVolume=0.0;                           //!< Ω
    ΔG_Map itsRho;                                  //!< the fit = ρ̃(G), private
    bool   itsFitted=false;
    projector_t itsProjector;                       //!< the projection half (ρ̃ on the same basis's overlap face)
};

} //namespace
