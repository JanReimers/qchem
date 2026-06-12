// File: BasisSet1/Atom/Evaluators/Internal/RKBL_Angular.C
// Mixin for RKB large-component angular ERIs. Holds κ and optional mj occupations
// and overrides CoulombAk/ExchangeAk via RelAngularIntegrals.
module;
#include <sstream>
export module qchem.BasisSet.Atom.Evaluators.Internal.RKBL_Angular;
export import qchem.BasisSet.Atom.Evaluators.IBS;

export namespace BasisSet::Atom::Evaluators
{

class RKBL_Angular : public virtual Evaluator
{
public:
    // mjs empty → full mj sum (closed shell or spin-averaged)
    RKBL_Angular(int _κ, const rvec_t& _mjs={}) : Evaluator(0), κ(_κ), mjs(_mjs) {}

    virtual rvec11_t    CoulombAk (const Evaluator& other) const override;
    virtual rvec11_t    ExchangeAk(const Evaluator& other) const override;
    virtual std::string AngularID () const override;

    int           Getκ  () const { return κ; }
    const rvec_t& Getmjs() const { return mjs; }

private:
    int    κ;
    rvec_t mjs;
};

} //namespace
