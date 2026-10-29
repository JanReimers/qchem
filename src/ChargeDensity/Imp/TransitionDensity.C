// File: ChargeDensity/Imp/TransitionDensity.C  The factory of the AO-matrix transition density.
module;
#include <memory>
#include <stdexcept>
#include <vector>
module qchem.ChargeDensity.TransitionDensity;
import qchem.ChargeDensity.Internal.TransitionDensity;   // the concrete

namespace qchem::ChargeDensity
{

template <class T> std::unique_ptr<TransitionDensity<T>>
AO_TransitionDensity_Factory(std::vector<TransitionBlock<T>> blocks, std::shared_ptr<const Symmetry::SelectionRule> rule)
{
    return std::make_unique<AO_TransitionDensity<T>>(std::move(blocks), std::move(rule));
}

namespace {
//! RestrictCoupling's view: every question but Coupling() is the inner density's.  The faces are cross-cast on the INNER
//! object (abstract -> abstract), so the view serves whatever the concrete serves and nothing more.
class CouplingView
    : public virtual TransitionDensity<dcmplx>
    , public virtual TransitionFourierDensity
    , public virtual ProjectableTransition
{
public:
    CouplingView(std::unique_ptr<TransitionDensity<dcmplx>> inner, std::shared_ptr<const Symmetry::SelectionRule> rule)
        : itsInner(std::move(inner)), itsRule(std::move(rule))
    {
        if (!itsInner || !itsRule) throw std::invalid_argument("RestrictCoupling: a null density or rule");
        const rvec3_t a=WaveVectorOf(*itsRule), b=WaveVectorOf(*itsInner->Coupling());
        if (a.x!=b.x || a.y!=b.y || a.z!=b.z)
            throw std::invalid_argument("RestrictCoupling: the restricted rule must keep the density's wave vector");
    }
    virtual std::shared_ptr<const Symmetry::SelectionRule> Coupling() const override {return itsRule;}
    virtual size_t Version() const override {return itsInner->Version();}   // the same content: a term's cache stays valid
    //! The channels are the INNER density's views: a term contracts a channel's CONTENT; only the kernel's pair loop, at the
    //! top level, reads the coupling.
    virtual const TransitionDensity<dcmplx>* Channel(const Spin& s) const override
    {return s==Spin::None ? this : itsInner->Channel(s);}
    virtual ΔGq_Map GetTransitionRepulsion(const BasisSet::cFIT_CD_ABS& fit) const override
    {
        auto* f=dynamic_cast<const TransitionFourierDensity*>(itsInner.get());
        if (!f) throw std::logic_error("RestrictCoupling: the inner density has no Hartree face");
        return f->GetTransitionRepulsion(fit);
    }
    virtual cvec_t ProjectOnto(const TransitionProjector& P) const override
    {
        auto* p=dynamic_cast<const ProjectableTransition*>(itsInner.get());
        if (!p) throw std::logic_error("RestrictCoupling: the inner density has no pair-projection face");
        return p->ProjectOnto(P);
    }
private:
    std::unique_ptr<TransitionDensity<dcmplx>>     itsInner;
    std::shared_ptr<const Symmetry::SelectionRule> itsRule;
};
} // namespace

std::unique_ptr<TransitionDensity<dcmplx>> RestrictCoupling(std::unique_ptr<TransitionDensity<dcmplx>> inner,
                                                            std::shared_ptr<const Symmetry::SelectionRule> rule)
{
    return std::make_unique<CouplingView>(std::move(inner), std::move(rule));
}

template std::unique_ptr<TransitionDensity<double>> AO_TransitionDensity_Factory<double>(std::vector<TransitionBlock<double>>, std::shared_ptr<const Symmetry::SelectionRule>);
template std::unique_ptr<TransitionDensity<dcmplx>> AO_TransitionDensity_Factory<dcmplx>(std::vector<TransitionBlock<dcmplx>>, std::shared_ptr<const Symmetry::SelectionRule>);

} // namespace
