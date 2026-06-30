module;
#include <cassert>
#include <iostream>
module qchem.BasisSet.Internal.ERI3;
import qchem.Blaze;

namespace qchem {

template <class T> double fnorm(const ERI3<T>& a, const ERI3<T>& b)
{
    double ret=0.0;
    assert(a.size()==b.size());
    auto mb=b.begin();
    for (auto ma:a)
    {
        double norm_ab=blazem::norm(ma-*mb);
        ret+=norm_ab*norm_ab;
        mb++;
    }
    return sqrt(ret);    
}

template <class T> double relative_fnorm(const ERI3<T>& a, const ERI3<T>& b)
{
    double ret=0.0;
    assert(a.size()==b.size());
    auto mb=b.begin();
    for (auto ma:a)
    {
        // std::cout << "ma=" << ma << std::endl << "mb=" << *mb << std::endl << std::endl ;
        double norm_ab=blazem::norm(ma-*mb);
        double avg_norm_ab=(blazem::norm(ma)+blazem::norm(*mb))/2.0;
        if (avg_norm_ab>0.0) norm_ab/=avg_norm_ab;
        ret+=norm_ab*norm_ab;
        mb++;
    }
    // std::cout << "--------------------------------------------------------------------------" << std::endl;
    return sqrt(ret);    
}

template double          fnorm<double>(const ERI3<double>& a, const ERI3<double>& b);
template double relative_fnorm<double>(const ERI3<double>& a, const ERI3<double>& b);

} // namespace qchem