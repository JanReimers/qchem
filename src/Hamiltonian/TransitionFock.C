// File: Hamiltonian/TransitionFock.C  The first-order Fock change of a linear response
// (doc/LinearResponsePlan.md §3 row C2, §3c).
//
// A DISTINCT TYPE FROM THE TRANSITION DENSITY (pin 20, ask what a matrix MEANS).  δD is contravariant and δF
// is covariant; their contraction is an energy.  One container for both would let a solver add a density to
// a potential -- as separate types that is a build error.
//
// ONE FACE ONLY, Matrix(bra, ket) -- the form a sum-over-states R0 consumes.  The Sternheimer form
// Apply(bra, ket, ψ) = δV|ψ> arrives as a cross-cast capability WITH the first Sternheimer concrete (ruling Q1,
// 2026-09-28: a face clause with no honest implementor is what R2.7 deleted from FittedCD).
module;
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <utility>
export module qchem.Hamiltonian.TransitionFock;
export import qchem.Symmetry.SelectionRule;
export import qchem.Symmetry.Irrep;
import qchem.Streamable;   // Irrep's operator<<
import qchem.Types;
import qchem.Blaze;

export namespace qchem::Hamiltonian
{

template <class T> class TransitionFock
{
public:
    virtual ~TransitionFock() = default;
    //! Which bra block each ket block couples to -- the SAME label as the transition density that induced it.
    virtual std::shared_ptr<const Symmetry::SelectionRule> Coupling() const = 0;
    //! δF on the coupled block pair (bra <- ket), in the AO basis: bra rows x ket columns.  THROWS on a pair
    //! this Fock change does not hold.
    virtual mat_t<T> Matrix(const Irrep& bra, const Irrep& ket) const = 0;
};

//! \brief δF as AO matrices, one per coupled block pair.  A VALUE: the response kernel fills it term by term
//! (Add accumulates), and so can a probe (a perturbation operator h^(1)).
template <class T> class AO_TransitionFock : public virtual TransitionFock<T>
{
    using rule_t=std::shared_ptr<const Symmetry::SelectionRule>;
public:
    explicit AO_TransitionFock(rule_t rule) : itsRule(std::move(rule))
    {
        if (!itsRule) throw std::invalid_argument("AO_TransitionFock: no selection rule");
    }
    virtual rule_t   Coupling() const override {return itsRule;}
    virtual mat_t<T> Matrix(const Irrep& bra, const Irrep& ket) const override
    {
        auto i=itsPairs.find({bra,ket});
        if (i==itsPairs.end())
        {
            std::ostringstream os;
            os << "AO_TransitionFock::Matrix: no block pair " << bra << " <- " << ket;
            throw std::out_of_range(os.str());
        }
        return i->second;
    }
    //! δF(bra <- ket) += \a m.  THROWS if the selection rule does not couple the pair, or the shape changes.
    void Add(const Irrep& bra, const Irrep& ket, const mat_t<T>& m)
    {
        if (!itsRule->Couples(*bra.sym, *ket.sym))
        {
            std::ostringstream os;
            os << "AO_TransitionFock::Add: the selection rule does not couple " << bra << " <- " << ket;
            throw std::invalid_argument(os.str());
        }
        auto [i,fresh]=itsPairs.try_emplace({bra,ket}, m);
        if (fresh) return;
        if (i->second.rows()!=m.rows() || i->second.columns()!=m.columns())
            throw std::invalid_argument("AO_TransitionFock::Add: a block pair's matrix changed shape");
        i->second+=m;
    }
private:
    //! Ordered by (bra, ket) -- Irrep is ordered by its sequence index, not a pointer.
    struct Less
    {
        bool operator()(const std::pair<Irrep,Irrep>& a, const std::pair<Irrep,Irrep>& b) const
        {return a.first<b.first || (!(b.first<a.first) && a.second<b.second);}
    };
    rule_t                                        itsRule;
    std::map<std::pair<Irrep,Irrep>,mat_t<T>,Less> itsPairs;
};

using rTransitionFock=TransitionFock<double>;   using cTransitionFock=TransitionFock<dcmplx>;

} // namespace
