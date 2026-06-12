// File: BasisSet1/Atom/Evaluators/Internal/NR_Angular.C
// Mixin for NR (non-relativistic) angular ERIs. Holds integer ml occupations
// and overrides CoulombAk/ExchangeAk via AngularIntegrals.
module;
#include <sstream>
export module qchem.BasisSet.Atom.Evaluators.Internal.NR_Angular;
export import qchem.BasisSet.Atom.Evaluators.IBS;

export namespace BasisSet::Atom::Evaluators
{

class NR_Angular : public virtual Evaluator
{
public:
    NR_Angular(const ivec_t& _mls) : Evaluator(0), mls(_mls) {}

    virtual rvec11_t    CoulombAk (const Evaluator& other) const override;
    virtual rvec11_t    ExchangeAk(const Evaluator& other) const override;
    virtual std::string AngularID () const override;

    const ivec_t& Getmls() const { return mls; }

private:
    ivec_t mls;
};

} //namespace
