// File: BasisSet/Atom/Evaluators/Internal/Imp/NR_Angular.C
module;
#include <blaze/Math.h>
#include <sstream>
module qchem.BasisSet.Atom.Evaluators.Internal.NR_Angular;
import qchem.BasisSet.Atom.Evaluators.Internal.AngularIntegrals;

namespace BasisSet::Atom::Evaluators
{

rvec11_t NR_Angular::DirectAk(const Evaluator& other) const
{
    const NR_Angular& o = dynamic_cast<const NR_Angular&>(other);
    int la=Getl(), lc=o.Getl();
    const ivec_t& amls=mls, cmls=o.mls;
    size_t nac=amls.size()*cmls.size();
    size_t g=size_t(2*la+1)*size_t(2*lc+1);
    rvec11_t Ak(0.0);
    if (nac==0 || nac==g)
    {
        Ak += AngularIntegrals::Direct (la,lc);
    }
    else
    {
        for (auto ma:amls)
        for (auto mc:cmls)
            Ak += AngularIntegrals::Direct (la,lc,ma,mc);
        Ak /= (double)nac;
    }
    return Ak;
}

rvec11_t NR_Angular::ExchangeAk(const Evaluator& other) const
{
    const NR_Angular& o = dynamic_cast<const NR_Angular&>(other);
    int la=Getl(), lb=o.Getl();
    const ivec_t& amls=mls, bmls=o.mls;
    size_t nab=amls.size()*bmls.size();
    size_t g=size_t(2*la+1)*size_t(2*lb+1);
    rvec11_t Ak(0.0);
    if (nab==0 || nab==g)
    {
        Ak += AngularIntegrals::Exchange(la,lb);
    }
    else
    {
        for (auto ma:amls)
        for (auto mb:bmls)
            Ak += AngularIntegrals::Exchange(la,lb,ma,mb);
        Ak /= (double)nab;
    }
    return Ak;
}

std::string NR_Angular::AngularID() const
{
    std::ostringstream os;
    os << Getl() << " {";
    for (auto ml:mls) os << ml << " ";
    os << "}";
    return os.str();
}

} //namespace
