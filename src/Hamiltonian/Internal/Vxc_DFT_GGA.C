// File: Hamiltonian/Internal/Vxc_DFT_GGA.C  The GGA exchange-correlation term (GGA stage 2, the molecular point route; GGA design note §4).
//
// E_xc = Σ_g w_g e(ρ↑, ρ↓, σ↑↑, σ↑↓, σ↓↓) on the Hamiltonian's nonlinear integration mesh, and its KS matrix per spin
//     H^σ_ij = Σ_g w_g [ v_ρσ χ_iχ_j + w_σ · ∇(χ_iχ_j) ],   w_σ = 2 v_σσσ ∇ρ_σ + v_σ↑↓ ∇ρ_σ̄,
// the exact derivative of that E (the sampler's {v, w} gather is the transpose of its (ρ, ∇ρ) sampling).
//
// TWO TERMS, TWO CONSTRUCTORS (user, 2026-11-20: the gradient is a TYPE).  The LDA term is tVxc_DFT, untouched.  This one DEMANDS a
// GGAFunctional and a density-gradient sampler -- concrete faces, no default, no null -- so its density state is {ρ↑, ρ↓, ∇ρ↑, ∇ρ↓}
// always and its field {v, w} always: there is no "gradients present only when" branch.  An unpolarized run is the ρ/2, ∇ρ/2 collapse the
// sampler applies (pin 5), fed to the same spin-native libxc call.
//
// Not yet here (named refusals elsewhere): forces and the virial (stage 5: the pullback needs ∇∇φ), the response kernel (stage 6).
// The Bloch axis (stage 3) is the same term on <dcmplx>, with the real TRIM blocks' cache beside it.
module;
#include <iosfwd>
#include <memory>
#include <optional>
#include <type_traits>
#include <complex>
#include <vector>
export module qchem.Hamiltonian.Internal.Vxc_DFT_GGA;
import qchem.Hamiltonian.Internal.Term;            // tDynamic_HT_Imp: the refresh-built block cache
import qchem.Hamiltonian.Internal.GGAFunctional;   // GGAFunctional, GGAInput / GGAOutput, RequireGradientSampled
import qchem.ChargeDensity.DensitySampler;         // tDensityGradientSampler: the (ρ, ∇ρ) quadrature and its {v, w} gather
import qchem.Hamiltonian.Types;
import qchem.Parallel;                             // LazyLock: the Fock blocks record the pin-28 source concurrently
import qchem.Types;

export namespace qchem::Hamiltonian
{

//! \brief WHAT ONLY A BLOCH RUN HAS, on \c tVxc_DFT_GGA<dcmplx> (GGA stage 3): the real TRIM block's cache (Step 3c), its blocks the periodic
//! sampler's real {v, w} gather.  No response face: a GGA response is stage 6, refused by name above.  (The LDA term's twin is
//! \c Vxc_DFT_BlochFaces; a virtual cannot carry a \c requires clause, so the body is the term's, constrained, and this routes to it.)
template <class D> class GGA_BlochFaces
    : public Dynamic_HT_RealBlock_Imp
{
protected:
    virtual rsmat_t MakeMatrix(const robs_t* bs, const Spin& s, DensityID id) const override
    {return static_cast<const D&>(*this).RealBlock(bs, s, id);}
    void BuildRealBlocks(const cbs_t* bs, const cChargeDensity* cd) const {Dynamic_HT_RealBlock_Imp::BuildBlocks(bs, cd);}
};
//! The molecular instantiation's stand-in: one block cache.
struct NoGGABlochFaces { template <class B, class C> void BuildRealBlocks(const B*, const C*) const {} };
template <class T, class D> using GGA_BlochBase = std::conditional_t<std::is_same_v<T,dcmplx>, GGA_BlochFaces<D>, NoGGABlochFaces>;

template <class T> class tVxc_DFT_GGA
    : public virtual tDynamic_HT<T>
    , protected      tDynamic_HT_Imp<T>
    , public         GGA_BlochBase<T,tVxc_DFT_GGA<T>>   // on <dcmplx>: the real TRIM block's cache
{
public:
    typedef std::shared_ptr<const GGAFunctional>                              gga_t;
    typedef std::shared_ptr<const ChargeDensity::tDensityGradientSampler<T>>  sampler_t;
    //! The real (molecular) term.  Both pointers are REQUIRED (THROWS \c invalid_argument on null).
    tVxc_DFT_GGA(gga_t xc, sampler_t sampler, SpinGroup g) requires std::is_same_v<T,double>;
    //! The Bloch (periodic) term: the same sampler, through its periodic face for the real TRIM blocks (GGA stage 3).
    tVxc_DFT_GGA(gga_t xc, std::shared_ptr<const ChargeDensity::DensityGradientSampler> sampler, SpinGroup g) requires std::is_same_v<T,dcmplx>;

    //! Warms ρ and ∇ρ for \a cd, and runs the forgotten-gradient detector on them (\c RequireGradientSampled): once per refresh, live.
    virtual void          RefreshForDensity(const tChargeDensity<T>* cd) const override;
    //! Builds the field {v, w} per channel for \a cd, then the blocks (which only read it).
    virtual void          BuildBlocks(const tbs_t<T>* bs, const tChargeDensity<T>* cd) const override;
    virtual void          GetEnergy(EnergyBreakdown&, const tDM_CD<T>*) const override;
    virtual std::vector<Spin> BlockSpins() const override {return SpinIrreps(itsGroup);}
    virtual std::ostream& Write(std::ostream&) const override;
protected:
    virtual hmat_t<T> MakeMatrix(const tobs_t<T>*, const Spin&, DensityID) const override;
    friend class GGA_BlochFaces<tVxc_DFT_GGA<T>>;
    //! A real TRIM block's {v, w} gather (Step 3c), through the periodic face.
    rsmat_t RealBlock(const robs_t* bs, const Spin& s, DensityID id) const requires std::is_same_v<T,dcmplx>;
private:
    //! The functional's input at my points for \a cd: the sampler's channels and the σ triple built from their gradients.
    GGAInput InputAt(const tChargeDensity<T>* cd) const;

    gga_t     itsXc;
    sampler_t itsSampler;   //!< the Hamiltonian's nonlinear integration mesh (V-NONLINEAR-MESH): I am one client of it
    std::shared_ptr<const ChargeDensity::DensityGradientSampler> itsPeriodic;   //!< the same sampler's periodic face (<dcmplx> only)
    SpinGroup itsGroup;
    //! What the Fock refresh built for one density, which every block then reads (V1.41 stage 2).
    struct Refreshed
    {
        DensityID  id=DensityID(static_cast<const tChargeDensity<T>*>(nullptr));
        rvec_t     v[2];
        rvec3vec_t w[2];
        std::optional<ChargeDensity::XCRhoSource> source;
        bool densityIsRhoD=false;
    };
    mutable Refreshed itsRefreshed;
    //! The pin-28 certificate of the last Fock build (as \c tVxc_DFT records it): what the energy pass reports.
    mutable std::optional<ChargeDensity::XCRhoSource> itsFockSource;
    mutable bool itsFockDensityIsRhoD=false;
    qchem::LazyLock itsFockSourceLock;   //!< every Fock block writes the source (the same value): concurrent blocks take turns
    //! The block path's read of \c itsRefreshed: throws if the refresh never saw the density; records the Fock pass's pin-28 source.
    const Refreshed& RefreshedFor(DensityID id) const;
};

} //namespace
