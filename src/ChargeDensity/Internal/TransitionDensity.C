// File: ChargeDensity/Internal/TransitionDensity.C  The AO-matrix transition density (q = 0) -- the concrete.
//
// INTERNAL: production code reaches it only through AO_TransitionDensity_Factory and the TransitionDensity
// face.  It is a module of its own (not hidden in the Imp unit) for ONE client: the finite-difference kernel
// oracle in src/Response/tests (ruling D6), which must read δD's AO blocks and is granted that through a
// friend declared in src/forward.H -- the face itself never hands the matrices out.
//
// COMPOSITION, NOT INHERITANCE (doc/LinearResponsePlan.md §3c).  It OWNS an ordinary tComposite_CD of δD
// leaves and forwards the one face a linear operator needs -- the HF sweep -- to it.  So J[δD] and K[δD] run
// through exactly the ground-state sweep code (Composite_HFSystem -> tHF_Pair_CD -> Orbital_HF_IBS), and no
// tDM_CD face is ever handed out.
//
// ⚠ THE LEAVES TAKE RhoRoute::Direct.  The factory's default route is the pivoted-Cholesky factor, which
// assumes D positive semi-definite; δD is traceless and INDEFINITE.  (Harmless for J/K, which never evaluate
// ρ(r), but R2's mesh sampling would inherit a factorisation that is invalid for it.)
module;
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include <vector>
#include "forward.H"   // TransitionDensityTests -- the FD-kernel oracle's friend (ruling D6)
export module qchem.ChargeDensity.Internal.TransitionDensity;
export import qchem.ChargeDensity.TransitionDensity;
import qchem.ChargeDensity;
import qchem.ChargeDensity.Factory;   // IrrepCD_Factory, RhoRoute
import qchem.CompositeCD;
import qchem.ChargeDensity.FourierDensity;   // the periodic face δρ̃ forwards (Hartree, R2)
import qchem.Fitting.FunctionFitter;          // Fitting::ScalarProjector (ProjectOnto)
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

//! The G-space face, periodic path only (Hartree reads δρ̃ through it, R2): the same forwarding.
template <class Self> class Transition_Fourier : public virtual FourierDensity
{
public:
    virtual ΔG_Map GetFourierDensity (const BasisSet::cFIT_SF_ABS& c) const override {return F().GetFourierDensity(c);}
    virtual ΔG_Map GetRepulsion3C    (const BasisSet::cFIT_CD_ABS& c) const override {return F().GetRepulsion3C(c);}
    virtual ΔG_Map GetRepulsion3C_Raw(const BasisSet::cFIT_CD_ABS& c) const override {return F().GetRepulsion3C_Raw(c);}
    virtual void   StarAverage(ΔG_Map& rg) const override {F().StarAverage(rg);}
    virtual rvec_t GetRhoOnGrid(const BasisSet::cFIT_SF_ABS& c) const override {return F().GetRhoOnGrid(c);}
private:
    const FourierDensity& F() const
    {
        auto* f=dynamic_cast<const FourierDensity*>(static_cast<const Self&>(*this).Operand());
        if (!f) throw std::logic_error("AO_TransitionDensity: the wrapped composite has no FourierDensity face");
        return *f;
    }
};
template <class T, class Self> using TransitionFourierBase =
    std::conditional_t<std::is_same_v<T,dcmplx>, Transition_Fourier<Self>, NoFourierDensity>;
template <class T, class Self> using TransitionHFBase =
    std::conditional_t<std::is_same_v<T,double>, Transition_HFSystem<Self>, NoTransition_HF>;

template <class T> class AO_TransitionDensity
    : public virtual TransitionDensity<T>
    , public virtual tProjectable_CD<T>                        //!< δρ sampled like any density (XC mesh, +U forward), R2
    , public TransitionHFBase<T, AO_TransitionDensity<T>>      //!< J/K: real path, R1
    , public TransitionFourierBase<T, AO_TransitionDensity<T>> //!< δρ̃ for Hartree: periodic path, R2
{
    using rule_t=std::shared_ptr<const Symmetry::SelectionRule>;
public:
    AO_TransitionDensity(std::vector<TransitionBlock<T>> blocks, rule_t rule)
        : itsOwned(std::make_unique<tComposite_CD<T>>()), itsCD(nullptr), itsRule(std::move(rule)), itsSpin(Spin::None)
    {
        if (!itsRule) throw std::invalid_argument("AO_TransitionDensity: no selection rule");
        if (blocks.empty()) throw std::invalid_argument("AO_TransitionDensity: no blocks");
        for (auto& b : blocks)
        {
            if (!b.bs) throw std::invalid_argument("AO_TransitionDensity: a block with no basis");
            if (!itsRule->Couples(*b.irrep.sym, *b.irrep.sym))
            {
                std::ostringstream os;
                os << "AO_TransitionDensity: block " << b.irrep << " is not coupled to ITSELF by this selection "
                      "rule -- a (k+q, k) pair transition density is stage R3's (doc/LinearResponsePlan.md §5)";
                throw std::invalid_argument(os.str());
            }
            itsOwned->Insert(std::unique_ptr<tDM_CD<T>>(IrrepCD_Factory<T>(b.dD, b.bs, b.irrep, RhoRoute::Direct)), b.irrep);
        }
        itsBlocks=std::move(blocks);
        itsCD=itsOwned.get();
        // Eager, so a const reader never races on a lazy build (the CompositeCD channel-view rule).
        for (Spin s : {Spin::Up, Spin::Down})
            if (const tChargeDensity<T>* ch=itsOwned->GetChannel(s))
                itsChannels[s]=std::unique_ptr<AO_TransitionDensity>(new AO_TransitionDensity(ch, itsRule, s));
    }

    virtual rule_t Coupling() const override {return itsRule;}
    virtual size_t Version () const override {return itsCD->Version();}
    virtual const TransitionDensity<T>* Channel(const Spin& s) const override
    {
        if (s==Spin::None) return this;
        if (itsOwned)                            // the whole density: its channel views
        {
            auto i=itsChannels.find(s);
            return i==itsChannels.end() ? nullptr : i->second.get();
        }
        return s==itsSpin ? this : nullptr;      // a view IS its own channel
    }
    //! \copydoc tProjectable_CD::ProjectOnto
    //! δρ projected exactly as the ground-state density is -- the leaves' own tables, linear in δD.
    virtual rvec_t ProjectOnto(const Fitting::ScalarProjector& p) const override
    {
        auto* pr=dynamic_cast<const tProjectable_CD<T>*>(itsCD);
        if (!pr) throw std::logic_error("AO_TransitionDensity: the wrapped composite cannot project");
        return pr->ProjectOnto(p);
    }
    //! What the capability faces forward to (Transition_HFSystem, Transition_Fourier).
    const tChargeDensity<T>* Operand() const {return itsCD;}

private:
    //! A channel VIEW: non-owning, over the parent composite's \a s channel.
    AO_TransitionDensity(const tChargeDensity<T>* channel, rule_t rule, Spin s)
        : itsCD(channel), itsRule(std::move(rule)), itsSpin(s) {}

    friend class ::TransitionDensityTests;        //!< the FD oracle reads itsBlocks (src/forward.H, D6)
    std::vector<TransitionBlock<T>>   itsBlocks;  //!< THE representation: δD per block (empty in a view)
    std::unique_ptr<tComposite_CD<T>> itsOwned;   //!< the same δD as IrrepCD leaves -- what the HF sweep runs on; null in a view
    const tChargeDensity<T>*          itsCD;      //!< the composite, or the channel view of it
    rule_t                            itsRule;
    Spin                              itsSpin;    //!< a view's channel (None for the whole density)
    std::map<Spin,std::unique_ptr<AO_TransitionDensity>> itsChannels;   //!< the Up/Down views (whole density only)
};

} // namespace
