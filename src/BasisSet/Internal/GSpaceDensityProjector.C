// File: BasisSet/Internal/GSpaceDensityProjector.C  The ρ̃ projector's CLASS (the face is BasisSet::PairProjector<dcmplx>).
//
// Internal by design (user 2026-11-08): every client only asks density.Project(p), a PairProjector<dcmplx> call, so only the
// factory (qchem.BasisSet.GSpaceDensityProjector) and this unit name the concrete type.  What it does: see the factory's header.
module;
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>
export module qchem.BasisSet.Internal.GSpaceDensityProjector;
export import qchem.BasisSet.PairProjector;     // PairProjector<dcmplx>, AnalyticFourier
import qchem.BasisSet.PairProjection;           // GSpacePairProjectionIntegrator
import qchem.BasisSet.Orbital_DFT_IBS;          // cFIT_SF_ABS, Orbital_DFT_IBS<U,dcmplx>
import qchem.BasisSet.G_FieldEvaluator;         // G_RasterTransform (the raster a field is sampled on)
import qchem.BasisSet.GMap;                     // GMapStarField (ρ̃ as a field the group averages)
import qchem.Symmetry;                          // SymMap
import qchem.Parallel;                          // LazyLock
import qchem.Reporting;                         // report::Timed

export namespace qchem::BasisSet
{

class GSpaceDensityProjector : public virtual PairProjector<dcmplx>
{
public:
    //! The constructor DEMANDS the fit basis \a sf (co-owned) and its raster face \a rt, of the same object (CLAUDE.md "Constructors
    //! pass concrete objects"); \c MakeGSpaceDensityProjector gathers them.  \a starAverageBucket names the ledger row the
    //! star-average is timed under (each client keeps its own, so the ledgers stay comparable across 3c).
    GSpaceDensityProjector(std::shared_ptr<const cFIT_SF_ABS> sf, const G_RasterTransform& rt, std::string starAverageBucket)
        : itsSF(std::move(sf)), itsRT(&rt), itsBucket(std::move(starAverageBucket))
    {
        if (!itsSF) throw std::logic_error("GSpaceDensityProjector: no fit basis");
        if (!itsSF->isOrtho())
            throw std::logic_error("GSpaceDensityProjector: the fit basis '"+itsSF->Name()+"' is not orthonormal, so its projection is not ρ̃");
    }

    const qcMesh::MatrixForward<double>& Forward(const Orbital_DFT_IBS<double,dcmplx>& orb) const override {return Integrator<double>(orb);}
    const qcMesh::MatrixForward<dcmplx>& Forward(const Orbital_DFT_IBS<dcmplx,dcmplx>& orb) const override {return Integrator<dcmplx>(orb);}
    //! A density with no D: its own Fourier components when it knows them, else sampled on my raster and transformed.
    Projections ProjectField(const ScalarFunction<double>& f) const override
    {
        // A cross-cast between abstract faces, and the one seam there is: a density with no D reaches me as a ScalarFunction
        // through tChargeDensity::Project's default, and a second Project entry per capability is what ruling (e) declined.
        if (auto* af=dynamic_cast<const AnalyticFourier*>(&f)) return Projections(af->FourierComponents(*itsSF));
        return Projections(itsRT->FieldCoeffs(itsRT->ForwardFFT(itsRT->Sample(f))));
    }
    size_t NumCoefficients() const override {return itsSF->GetNumFunctions();}
    //! A zero density's ρ̃: the empty map (an absent G reads as 0).
    Projections Zero() const override {return Projections(ΔG_Map{});}
    //! The IBZ star-average of ρ̃: my functions are plane waves, so a point op permutes their G's.
    void StarAverage(Projections& p, const Symmetry::StarGroup& group) const override
    {
        if (group.IsTrivial()) return;
        qchem::report::Timed timed(itsBucket);
        ΔG_Map rho=std::move(p).GMap();   // moved out, averaged in place, moved back
        GMapStarField f(rho);
        group.StarAverage(f);
        p=Projections(std::move(rho));
    }

private:
    //! The projection VIEW for one orbital block, over the basis's cached OVERLAP tensor, by spatial symmetry (as every fitter).
    template <class U> const GSpacePairProjectionIntegrator<U>& Integrator(const Orbital_DFT_IBS<U,dcmplx>& orb) const
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

    std::shared_ptr<const cFIT_SF_ABS> itsSF;     //!< the plane-wave fit basis (its overlap face)
    const G_RasterTransform*           itsRT;     //!< ...its raster
    std::string                        itsBucket; //!< the star-average's ledger row
    mutable SymMap<GSpacePairProjectionIntegrator<double>> itsIntR;   //!< real TRIM blocks
    mutable SymMap<GSpacePairProjectionIntegrator<dcmplx>> itsIntC;   //!< Bloch blocks
    qchem::LazyLock itsCacheLock;
};

} //namespace
