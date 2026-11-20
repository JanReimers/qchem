// File: BasisSet/Imp/RasterDensityProjector.C  The raw ρ_DM projector's class and factory (only this unit names the concrete type).
module;
#include <memory>
#include <mutex>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
module qchem.BasisSet.RasterDensityProjector;
import qchem.BasisSet.PairProjection;      // PairProjection<dcmplx>::applyRaw (the block's raw collocation)
import qchem.BasisSet.G_FieldEvaluator;    // G_RasterTransform (the raster's size)
import qchem.BasisSet.G_RasterStarField;   // the real-space star field over a raster array
import qchem.Symmetry;                     // StarGroup, StarField
import qchem.SlotMap;                      // InstanceSlotMap: the per-block forwards, by the block's instance
import qchem.Parallel;                     // LazyLock
import qchem.Mesh.Integrator;              // MatrixForward<U>

namespace qchem::BasisSet
{

namespace {

//! One block's RAW forward: \f$D\mapsto\rho_{DM}(r_g)\f$ on the raster, through the block's own collocation (\c applyRaw).  The
//! factored overload is the face's default (form \f$LL^\dagger\f$, delegate): a raster collocation has no thin-factor kernel.
template <class U> class RawRasterForward : public virtual qcMesh::MatrixForward<U>
{
public:
    RawRasterForward(const PairProjection<dcmplx>& g, size_t n) : itsG(&g), itsN(n) {}
    using qcMesh::MatrixForward<U>::Forward;   // keep the factored overload visible (qchem.Mesh.Integrator's warning)
    Projections Forward(const hmat_t<U>& D) const override {return Projections(itsG->applyRaw(D));}
    size_t NumCoefficients() const override {return itsN;}
private:
    const PairProjection<dcmplx>* itsG;   //!< borrowed: owned by the orbital basis's integral cache
    size_t itsN;                          //!< the raster's point count
};

class RasterDensityProjector : public virtual PairProjector<dcmplx>
{
public:
    RasterDensityProjector(std::shared_ptr<const cFIT_SF_ABS> sf, const G_RasterTransform& rt, const G_RasterStarField& star)
        : itsSF(std::move(sf)), itsRT(&rt), itsStar(&star) {}

    const qcMesh::MatrixForward<double>& Forward(const Orbital_DFT_IBS<double,dcmplx>& orb) const override {return Raw<double>(orb);}
    const qcMesh::MatrixForward<dcmplx>& Forward(const Orbital_DFT_IBS<dcmplx,dcmplx>& orb) const override {return Raw<dcmplx>(orb);}
    //! A density with no D: its field's own raster values (a mixer's shadow, or none for a seed).  A field without the face has no
    //! raster values to give, and sampling it here would quietly turn "raw" into "sampled" (pin 28): THROWS.
    Projections ProjectField(const ScalarFunction<double>& f) const override
    {
        // A cross-cast between abstract faces, the one seam there is (as the ρ̃ projector's AnalyticFourier).
        if (auto* rv=dynamic_cast<const RasterValuedField*>(&f)) return Projections(rv->RasterValues(*itsSF));
        throw std::logic_error("RasterDensityProjector: a density with no density matrix and no raster values of its own (RasterValuedField) "
                               "has no RAW rho_DM on the raster of '"+itsSF->Name()+"'");
    }
    size_t NumCoefficients() const override {return itsRT->RasterSize();}
    //! The IBZ star-average IN REAL SPACE: my raster says how its voxels transform under one op (a permutation plus an FFT glide), the
    //! composite's group iterates the ops.  Exactly the average the k-composite applied to its summed raster before 2026-11-19.
    void StarAverage(Projections& p, const Symmetry::StarGroup& group) const override
    {
        rvec_t rho=std::move(p).Vector();   // moved out, averaged in place, moved back
        {
            const std::unique_ptr<Symmetry::StarField> f=itsStar->MakeStarField(rho);
            group.StarAverage(*f);
        }
        p=Projections(std::move(rho));
    }

private:
    template <class U> const RawRasterForward<U>& Raw(const Orbital_DFT_IBS<U,dcmplx>& orb) const
    {
        auto& cache=Cache<U>();
        std::lock_guard<std::mutex> lock(itsCacheLock.m);
        return cache.GetOrFill(orb.Instance(), [&]
        {
            const PairProjection<dcmplx>& g=orb.Overlap3C(*itsSF);
            if (!g.applyRaw)
                throw std::logic_error("RasterDensityProjector: the orbital block '"+orb.Name()+"' has no raw collocation on the raster of '"
                                       +itsSF->Name()+"' -- the setup built a raster projector for a lineage without the raw route "
                                       "(RasterPairCollocation)");
            return RawRasterForward<U>(g, itsRT->RasterSize());
        });
    }
    template <class U> InstanceSlotMap<RawRasterForward<U>>& Cache() const
    { if constexpr (std::is_same_v<U,double>) return itsRawR; else return itsRawC; }

    std::shared_ptr<const cFIT_SF_ABS> itsSF;     //!< the plane-wave fit basis whose raster this is
    const G_RasterTransform*           itsRT;     //!< ...its raster
    const G_RasterStarField*           itsStar;   //!< ...and how the raster's voxels transform
    mutable InstanceSlotMap<RawRasterForward<double>> itsRawR;   //!< real TRIM blocks
    mutable InstanceSlotMap<RawRasterForward<dcmplx>> itsRawC;   //!< Bloch blocks
    qchem::LazyLock itsCacheLock;
};

} //anon

std::unique_ptr<PairProjector<dcmplx>> MakeRasterDensityProjector(std::shared_ptr<const cFIT_SF_ABS> sf)
{
    if (!sf) throw std::logic_error("MakeRasterDensityProjector: no fit basis");
    // The factory gathers the two faces of the same object (CLAUDE.md "constructors pass concrete objects").
    auto* rt=dynamic_cast<const G_RasterTransform*>(sf.get());
    if (!rt) throw std::logic_error("MakeRasterDensityProjector: the fit basis '"+sf->Name()+"' has no raster (G_RasterTransform)");
    auto* star=dynamic_cast<const G_RasterStarField*>(sf.get());
    if (!star) throw std::logic_error("MakeRasterDensityProjector: the fit basis '"+sf->Name()+"' cannot star-average its raster (G_RasterStarField)");
    return std::make_unique<RasterDensityProjector>(std::move(sf), *rt, *star);
}

} //namespace
