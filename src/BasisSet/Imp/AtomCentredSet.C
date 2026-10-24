// File: BasisSet/Imp/AtomCentredSet.C  The atom-centred set's value type: radial normalisation.
module;
#include <cmath>
#include <stdexcept>
#include <vector>
module qchem.BasisSet.AtomCentredSet;

namespace qchem::BasisSet
{

// int_0^inf r^{2l+2+2(n+n')} e^{-(a+a') r^2} dr = Gamma(p)/(2 (a+a')^p),  p = l + n + n' + 3/2.
double GaussianAtomCentredSet::RadialNorm2(int l, const std::vector<Math::Gaussian>& radial)
{
    double s=0.0;
    for (const auto& a : radial)
        for (const auto& b : radial)
        {
            const double p=l+a.n+b.n+1.5;
            s+=a.c*b.c*std::tgamma(p)/(2.0*std::pow(a.alpha+b.alpha,p));
        }
    return s;
}

size_t GaussianAtomCentredSet::Add(const rvec3_t& centre, int l, std::vector<Math::Gaussian> radial)
{
    if (l<0) throw std::invalid_argument("GaussianAtomCentredSet::Add: l must be >= 0");
    if (radial.empty()) throw std::invalid_argument("GaussianAtomCentredSet::Add: empty radial");
    const double n2=RadialNorm2(l,radial);
    if (!(n2>0.0)) throw std::invalid_argument("GaussianAtomCentredSet::Add: radial has zero norm");
    const double s=1.0/std::sqrt(n2);
    for (auto& t : radial) t.c*=s;
    itsM.push_back({centre,l,std::move(radial)});
    return itsM.size()-1;
}

} //namespace
