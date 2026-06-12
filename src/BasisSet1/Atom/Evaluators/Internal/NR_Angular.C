// File: BasisSet1/Atom/Evaluators/Internal/NR_Angular.C
// Mixin for NR (non-relativistic) angular ERIs. Holds integer ml occupations
// and overrides CoulombAk/ExchangeAk via AngularIntegrals.
module;
#include <sstream>
export module qchem.BasisSet.Atom.Evaluators.Internal.NR_Angular;
export import qchem.BasisSet.Atom.Evaluators.IBS;
import qchem.Symmetry.Spherical;

export namespace BasisSet::Atom::Evaluators
{

class NR_Angular : public virtual Evaluator
{
public:
    NR_Angular(int _l, const ivec_t& _mls) : Evaluator(_l), l(_l), mls(_mls) {}
    NR_Angular(const sym_t& sym) : NR_Angular(Symmetry::Getl(sym),Symmetry::Getmls(sym)) {};
    virtual rvec11_t    CoulombAk (const Evaluator& other) const override;
    virtual rvec11_t    ExchangeAk(const Evaluator& other) const override;
    virtual std::string AngularID () const override;

    virtual int           Getl  () const override {return l;}
    virtual const ivec_t& Getmls() const { return mls; }

private:
    int    l;
    ivec_t mls;
};

} //namespace
