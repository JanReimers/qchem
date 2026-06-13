// File: BasisSet1/Atom/Evaluators/Internal/Imp/RKBL_Angular.C
module;
#include <blaze/Math.h>
#include <sstream>
#include <cassert>
module qchem.BasisSet.Atom.Evaluators.Internal.RKBL_Angular;
import qchem.BasisSet.Atom.Evaluators.Internal.RelAngularIntegrals;
import qchem.Symmetry.Spherical; // SphericalSpinor::j(κ)

namespace BasisSet::Atom::Evaluators
{

rvec11_t RKB_Angular::CoulombAk(const Evaluator& other) const
{
    const RKB_Angular& o = dynamic_cast<const RKB_Angular&>(other);
    double ja = Symmetry::SphericalSpinor::j(κ);
    double jc = Symmetry::SphericalSpinor::j(o.κ);
    size_t ga = (size_t)(2*ja+1);
    size_t gc = (size_t)(2*jc+1);
    size_t nac = mjs.size() * o.mjs.size();
    rvec11_t Ak(0.0);
    if (nac==0 || nac==ga*gc)
    {
        // Full mj sum; divide by degeneracy to match the per-pair NR convention.
        Ak += RelAngularIntegrals::Coulomb(κ, o.κ);
        Ak /= (double)(ga*gc);
    }
    else
    {
        for (double mja : mjs)
        for (double mjc : o.mjs)
            Ak += RelAngularIntegrals::Coulomb(κ, o.κ, mja, mjc);
        Ak /= (double)nac;
    }
    return Ak;
}

rvec11_t RKB_Angular::ExchangeAk(const Evaluator& other) const
{
    const RKB_Angular& o = dynamic_cast<const RKB_Angular&>(other);
    double ja = Symmetry::SphericalSpinor::j(κ);
    double jb = Symmetry::SphericalSpinor::j(o.κ);
    size_t ga = (size_t)(2*ja+1);
    size_t gb = (size_t)(2*jb+1);
    size_t nab = mjs.size() * o.mjs.size();
    rvec11_t Ak(0.0);
    if (nab==0 || nab==ga*gb)
    {
        Ak += RelAngularIntegrals::Exchange(κ, o.κ);
        Ak /= (double)(ga*gb);
    }
    else
    {
        for (double mja : mjs)
        for (double mjb : o.mjs)
            Ak += RelAngularIntegrals::Exchange(κ, o.κ, mja, mjb);
        Ak /= (double)nab;
    }
    return Ak;
}

std::string RKB_Angular::AngularID() const
{
    std::ostringstream os;
    os << "κ=" << κ << " {";
    for (auto mj:mjs) os << mj << " ";
    os << "}";
    return os.str();
}

} //namespace
