// File: BasisSet/Radial/Evaluators/Internal/NR_Angular.C
// Mixin for NR (non-relativistic) angular ERIs. Holds integer ml occupations
// and overrides DirectAk/ExchangeAk via AngularIntegrals.
module;
#include <sstream>
export module qchem.BasisSet.Radial.Evaluators.Internal.NR_Angular;
export import qchem.BasisSet.Radial.Evaluators;
import qchem.Symmetry.Atom.Spherical;
import qchem.BasisSet.Radial.Evaluators.Internal.Wigner3j;

export namespace qchem::BasisSet::Radial::Evaluators
{

class NR_Angular : public virtual Angular
{
public:
    //! Throws std::out_of_range for l beyond the 3j table (Wigner3j::CheckL): checked once here, not per lookup.
    NR_Angular(int _l, const ivec_t& _mls) : l(_l), mls(_mls) {Wigner3j::CheckL(l);}
    NR_Angular(const sym_t& sym) : NR_Angular(Symmetry::Atom::Getl(sym),Symmetry::Atom::Getmls(sym)) {};
    virtual rvec11_t    DirectAk  (const Evaluator& other) const override;
    virtual rvec11_t    ExchangeAk(const Evaluator& other) const override;
    virtual std::string AngularID () const override;

    virtual int           Getl  () const {return l;}
    virtual const ivec_t& Getmls() const {return mls; }

private:
    int    l;
    ivec_t mls;
};

} //namespace
