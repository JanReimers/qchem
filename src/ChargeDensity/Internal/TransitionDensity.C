// File: ChargeDensity/Internal/TransitionDensity.C  The AO-matrix transition density -- the concrete.
//
// INTERNAL: production code reaches it only through AO_TransitionDensity_Factory and the TransitionDensity
// face.  It is a module of its own (not hidden in the Imp unit) for ONE client: the finite-difference kernel
// oracle in src/Response/tests (ruling D6), which must read δD's AO blocks and is granted that through a
// friend declared in src/forward.H -- the face itself never hands the matrices out.
//
// TWO PATHS, ONE REPRESENTATION (a list of (k+q, k) block pairs with δD on each; doc/LinearResponsePlan.md §3d C1):
//  * REAL (T = double, the molecular q = 0 path, R1): COMPOSITION, NOT INHERITANCE.  It OWNS an ordinary
//    tComposite_CD of δD leaves and forwards the one face a linear operator needs -- the HF sweep -- to it, so
//    J[δD] and K[δD] run through exactly the ground-state sweep code and no tDM_CD face is ever handed out.
//    ⚠ THE LEAVES TAKE RhoRoute::Direct: the factory's default route is the pivoted-Cholesky factor, which
//    assumes D positive semi-definite; δD is traceless and INDEFINITE.
//  * PERIODIC (T = dcmplx, any q): it holds the pairs and NOTHING ELSE, and answers its two faces by contracting
//    each pair through the KET block's B2 capabilities (Hartree) or the sampler's projector (XC).  ★ It no longer
//    forwards FourierDensity / tProjectable_CD through a composite of IrrepCD leaves (R2's route, retired at R3
//    step 3): those leaves STAR-AVERAGED δρ on an imposed run and folded δD through the T3 stream fold, and a
//    perturbation breaks the imposed group (§3d finding 5).  New code with no fold in it closes that by
//    construction, and q = 0 goes through it too ("q = 0 is the special case", §5).
module;
#include <algorithm>   // std::max (the HF-path symmetry check)
#include <cmath>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <stdexcept>
#include <type_traits>
#include <vector>
#include "forward.H"   // TransitionDensityTests -- the FD-kernel oracle's friend (ruling D6)
export module qchem.ChargeDensity.Internal.TransitionDensity;
export import qchem.ChargeDensity.TransitionDensity;
import qchem.ChargeDensity;
import qchem.ChargeDensity.Factory;   // IrrepCD_Factory, RhoRoute
import qchem.CompositeCD;
import qchem.BasisSet.Transition_DFT_IBS;     // the KET block's B2 pair capability (Hartree), cross-cast to
import qchem.Streamable;              // Irrep's operator<< (the failure message names the block)

export namespace qchem::ChargeDensity
{

//! The HF sweep face, real path only: forward to the composite (or its channel view) this density wraps.
//! CRTP for the same reason as Composite_HFSystem: the complex instantiation must not even DECLARE members
//! for an operation the periodic path never has.
template <class Self> class Transition_HFSystem : public virtual tHF_System_CD<double>
{
public:
    virtual void AccumulateDirectAll  (std::vector<hmat_t<double>>& Jall) const override {Sweep().AccumulateDirectAll(Jall);}
    virtual void AccumulateExchangeAll(std::vector<hmat_t<double>>& Kall) const override {Sweep().AccumulateExchangeAll(Kall);}
private:
    const tHF_System_CD<double>& Sweep() const
    {
        // abstract -> abstract: a real composite (or a channel view of one) carries the sweep face.
        auto* s=dynamic_cast<const tHF_System_CD<double>*>(static_cast<const Self&>(*this).Operand());
        if (!s) throw std::logic_error("AO_TransitionDensity: the wrapped composite has no HF sweep face");
        return *s;
    }
};
struct NoTransition_HF {};

//! The two PERIODIC faces (T = dcmplx), CRTP over the concrete's pairs and wave vector: the Hartree field summed
//! over the pairs through each KET block's B2 capability, and the sampler projection summed the same way.
template <class Self> class Transition_Periodic
    : public virtual TransitionFourierDensity
    , public virtual ProjectableTransition
{
public:
    virtual ΔGq_Map GetTransitionRepulsion(const BasisSet::cFIT_CD_ABS& fit) const override
    {
        const Self& me=static_cast<const Self&>(*this);
        ΔGq_Map out;
        out.q=me.WaveVector();
        for (const auto* p : me.Pairs())
        {
            const ΔGq_Map v=Ket(*p).TransitionRepulsion(fit, Bra(*p), p->dD, out.q);
            for (const auto& [dm,c] : v.c) out.c[dm]+=c;
        }
        return out;
    }
    virtual cvec_t ProjectOnto(const TransitionProjector& P) const override
    {
        const Self& me=static_cast<const Self&>(*this);
        cvec_t out(P.NumCoefficients(), dcmplx(0.0));
        for (const auto* p : me.Pairs()) out+=P.Forward(*p->braBs, *p->ketBs, p->dD, me.WaveVector());
        return out;
    }
private:
    // abstract -> abstract: a collocating Bloch block carries the pair face; a basis without it fails loudly.
    static const BasisSet::Transition_DFT_IBS& Face(const tobs_t<dcmplx>* bs, const char* which)
    {
        auto* t=dynamic_cast<const BasisSet::Transition_DFT_IBS*>(bs);
        if (!t) throw std::logic_error(std::string("AO_TransitionDensity: the ")+which+" block has no Transition_DFT_IBS "
                                       "face -- this basis cannot collocate a (k+q, k) transition density");
        return *t;
    }
    static const BasisSet::Transition_DFT_IBS& Ket(const TransitionBlock<dcmplx>& p) {return Face(p.ketBs, "ket");}
    static const BasisSet::Transition_DFT_IBS& Bra(const TransitionBlock<dcmplx>& p) {return Face(p.braBs, "bra");}
};
struct NoTransition_Periodic {};
template <class T, class Self> using TransitionPeriodicBase =
    std::conditional_t<std::is_same_v<T,dcmplx>, Transition_Periodic<Self>, NoTransition_Periodic>;
template <class T, class Self> using TransitionHFBase =
    std::conditional_t<std::is_same_v<T,double>, Transition_HFSystem<Self>, NoTransition_HF>;

template <class T> class AO_TransitionDensity
    : public virtual TransitionDensity<T>
    , public TransitionHFBase<T, AO_TransitionDensity<T>>        //!< J/K: real path, q = 0 (R1)
    , public TransitionPeriodicBase<T, AO_TransitionDensity<T>>  //!< δV_H(G+q) + the XC sampler projection: periodic path (R3)
{
    using rule_t=std::shared_ptr<const Symmetry::SelectionRule>;
    using pairs_t=std::vector<const TransitionBlock<T>*>;
public:
    AO_TransitionDensity(std::vector<TransitionBlock<T>> blocks, rule_t rule)
        : itsCD(nullptr), itsRule(std::move(rule)), itsSpin(Spin::None), itsVersion(NextDensityVersion())
    {
        if (!itsRule) throw std::invalid_argument("AO_TransitionDensity: no selection rule");
        if (blocks.empty()) throw std::invalid_argument("AO_TransitionDensity: no blocks");
        // The ONE wave vector of every pair (a rule without the face is q = 0 as far as a lattice is concerned).
        itsQ=WaveVectorOf(*itsRule);
        for (auto& b : blocks)
        {
            if (!b.braBs || !b.ketBs) throw std::invalid_argument("AO_TransitionDensity: a block pair with no basis");
            std::ostringstream os;
            if (!itsRule->Couples(*b.bra.sym, *b.ket.sym))
                os << "AO_TransitionDensity: this selection rule does not couple ket block " << b.ket << " to bra block " << b.bra;
            else if (b.bra.ms!=b.ket.ms)
                os << "AO_TransitionDensity: a spin-flip pair " << b.ket << " -> " << b.bra << " (not built)";
            else if (b.dD.rows()!=b.braBs->GetNumFunctions() || b.dD.columns()!=b.ketBs->GetNumFunctions())
                os << "AO_TransitionDensity: δD on " << b.ket << " -> " << b.bra << " is not bra x ket";
            if (!os.str().empty()) throw std::invalid_argument(os.str());
        }
        itsBlocks=std::move(blocks);
        for (const auto& b : itsBlocks) itsPairs.push_back(&b);   // stable: itsBlocks never changes after this
        if constexpr (std::is_same_v<T,double>)
        {
            // The HF sweep path: q = 0 leaves over a composite, exactly the ground-state sweep.
            itsOwned=std::make_unique<tComposite_CD<T>>();
            for (const auto& b : itsBlocks)
                itsOwned->Insert(std::unique_ptr<tDM_CD<T>>(IrrepCD_Factory<T>(HermitianBlock(b), b.ketBs, b.ket, RhoRoute::Direct)), b.ket);
            itsCD=itsOwned.get();
        }
        // Eager, so a const reader never races on a lazy build (the CompositeCD channel-view rule).
        for (Spin s : {Spin::Up, Spin::Down})
        {
            pairs_t sub;
            for (const auto* p : itsPairs) if (p->ket.ms==s) sub.push_back(p);
            if (sub.empty()) continue;
            const tChargeDensity<T>* ch=nullptr;
            if constexpr (std::is_same_v<T,double>)
                if (!(ch=itsOwned->GetChannel(s))) continue;
            itsChannels[s]=std::unique_ptr<AO_TransitionDensity>(new AO_TransitionDensity(ch, std::move(sub), itsRule, s, itsVersion, itsQ));
        }
    }

    virtual rule_t Coupling() const override {return itsRule;}
    virtual size_t Version () const override {return itsVersion;}
    virtual const TransitionDensity<T>* Channel(const Spin& s) const override
    {
        if (s==Spin::None) return this;
        if (itsSpin==Spin::None)                 // the whole density: its channel views
        {
            auto i=itsChannels.find(s);
            return i==itsChannels.end() ? nullptr : i->second.get();
        }
        return s==itsSpin ? this : nullptr;      // a view IS its own channel
    }
    //! What the HF sweep face forwards to (real path; null on the periodic path).
    const tChargeDensity<T>* Operand   () const {return itsCD;}
    //! The pairs of this density (a view: of its channel) and the one wave vector they share -- for the CRTP faces.
    const pairs_t&           Pairs     () const {return itsPairs;}
    const rvec3_t&           WaveVector() const {return itsQ;}

private:
    //! A channel VIEW: non-owning, over the parent's \a s pairs (and, on the real path, its composite's channel).
    AO_TransitionDensity(const tChargeDensity<T>* channel, pairs_t pairs, rule_t rule, Spin s, size_t version, rvec3_t q)
        : itsCD(channel), itsRule(std::move(rule)), itsSpin(s), itsVersion(version), itsQ(q), itsPairs(std::move(pairs)) {}
    //! The real path's leaf: q = 0 needs bra == ket and a Hermitian δD (R0 maps Hermitian to Hermitian, so
    //! anything else is a defect upstream -- THROWN, not silently symmetrised).
    static hmat_t<T> HermitianBlock(const TransitionBlock<T>& b)
    {
        if (b.bra.SequenceIndex()!=b.ket.SequenceIndex())
            throw std::invalid_argument("AO_TransitionDensity: a (bra != ket) pair on the REAL (HF sweep) path -- q != 0 "
                                        "HF response is not built");
        const size_t n=b.dD.rows();
        double scale=0.0, anti=0.0;
        for (size_t i=0;i<n;i++) for (size_t j=0;j<n;j++)
        {
            scale=std::max(scale, std::abs(b.dD(i,j)));
            anti =std::max(anti,  std::abs(b.dD(i,j)-b.dD(j,i)));
        }
        if (anti>1e-10*scale+1e-300)
            throw std::invalid_argument("AO_TransitionDensity: δD on the HF sweep path is not symmetric");
        hmat_t<T> H(n);
        for (size_t i=0;i<n;i++) for (size_t j=i;j<n;j++) H(i,j)=b.dD(i,j);
        return H;
    }

    friend class ::TransitionDensityTests;        //!< the FD oracle reads itsBlocks (src/forward.H, D6)
    std::vector<TransitionBlock<T>>   itsBlocks;  //!< THE representation: δD per pair (empty in a view)
    std::unique_ptr<tComposite_CD<T>> itsOwned;   //!< real path: the same δD as IrrepCD leaves (the HF sweep); null in a view
    const tChargeDensity<T>*          itsCD;      //!< real path: the composite, or the channel view of it
    rule_t                            itsRule;
    Spin                              itsSpin;    //!< a view's channel (None for the whole density)
    size_t                            itsVersion; //!< one serial for the density and its views (NextDensityVersion)
    rvec3_t                           itsQ{0,0,0};//!< the wave vector (fractional), from the rule
    pairs_t                           itsPairs;   //!< the pairs this density (or view) answers for
    std::map<Spin,std::unique_ptr<AO_TransitionDensity>> itsChannels;   //!< the Up/Down views (whole density only)
};

} // namespace
