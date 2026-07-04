// File: BasisSet/Atom/Evaluators/Internal/RKBL_Angular.C
// Mixin for RKB large-component angular ERIs. Holds κ and optional mj occupations
// and overrides DirectAk/ExchangeAk via RelAngularIntegrals.
module;
#include <sstream>
export module qchem.BasisSet.Atom.Evaluators.Internal.RKBL_Angular;
export import qchem.BasisSet.Atom.Evaluators;
import qchem.Symmetry.Atom.Spherical;

export namespace qchem::BasisSet::Atom::Evaluators
{

class RKB_Angular : public virtual Angular
{
public:
    // mjs empty → full mj sum (closed shell or spin-averaged)
    RKB_Angular(int _κ, const rvec_t& _mjs={}) : κ(_κ), mjs(_mjs) {}
    RKB_Angular(const sym_t& sym)
        : κ  (Symmetry::Getκ  (sym))
        , mjs(Symmetry::Getmjs(sym)) 
        {}

    virtual rvec11_t    DirectAk  (const Evaluator& other) const override;
    virtual rvec11_t    ExchangeAk(const Evaluator& other) const override;
    virtual std::string AngularID () const override;

    virtual int   Getl  () const { return Symmetry::Atom::SphericalSpinor::l(κ);}
    int           Getκ  () const { return κ; }
    const rvec_t& Getmjs() const { return mjs; }

private:
    int    κ;
    rvec_t mjs;
};

} //namespace
