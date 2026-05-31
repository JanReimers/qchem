// File:: BasisSet1/Atom/Evaluators/Imp/Evaluator.C
module;
#include <blaze/Math.h>

module qchem.BasisSet.Atom.Evaluators.IBS;
import qchem.BasisSet.Atom.Evaluators.Internal.AngularIntegrals;

namespace BasisSet::Atom::Evaluators
{
Evaluator::Evaluator(const Irrep_QNs::sym_t& y) :  l(0), mls({}),ns(0),grouper(0)
{
    const Yl_Sym* yl=dynamic_cast<const Yl_Sym*>(y.get());
    assert(yl);
    l=yl->GetL();
    const Ylm_Sym* ylm=dynamic_cast<const Ylm_Sym*>(y.get());
    if (ylm)
    {
        mls=ylm->Getmls();
    }
}
std::string Evaluator::AngularID() const
{
     std::ostringstream os;
     os << l << " {";
     for (auto ml:mls) os << ml << " ";
     os << "}";
     return os.str();
}

rvec11_t Evaluator::Coulomb_AngularIntegrals(const Evaluator& a,const Evaluator& c)
{
    rvec11_t Ak(0.0);
    int la=a.Getl(),lc=c.Getl();
    const Evaluator::ivec_t& amls=a.Getmls(),cmls=c.Getmls();
    size_t nac=amls.size()*cmls.size();
    size_t g=(2*la+1)*(2*lc+1); //degenracy
    if (nac==g || nac==0)
    {
        Ak+=AngularIntegrals::Coulomb(la,lc);
    }
    else
    {
        // For direct integrals these actually factor.  But for exchange they do not.
        // So it may not be worth while coding the factored version here.
        for (auto ma:amls)
        for (auto mc:cmls)
            Ak+=AngularIntegrals::Coulomb(la,lc,ma,mc);
        Ak/=(double)nac;
    }
    return Ak;
}
rvec11_t Evaluator::ExchangeAngularIntegrals(const Evaluator& a,const Evaluator& b)
{
    rvec11_t Ak(0.0);
    int la=a.Getl(),lb=b.Getl();
    const Evaluator::ivec_t& amls=a.Getmls(),bmls=b.Getmls();
    size_t nab=amls.size()*bmls.size();
    size_t g=(2*la+1)*(2*lb+1); //degenracy
    if (nab==0 || nab==g)
    {
        Ak+=AngularIntegrals::Exchange(la,lb);
    }
    else
    {
        for (auto ma:amls)
        for (auto mb:bmls)
            Ak+=AngularIntegrals::Exchange(la,lb,ma,mb);
        Ak/=nab;
    }
    return Ak;
}    

} //namespace