module;
#include <cassert>
#include <iostream>
module qchem.BasisSet.Internal.Projector3;
import qchem.Blaze;

namespace qchem {

template <class T> double fnorm(const Projector3<T>& a, const Projector3<T>& b)
{
    double ret=0.0;
    assert(a.dense.size()==b.dense.size());
    auto mb=b.dense.begin();
    for (auto ma:a.dense)
    {
        double norm_ab=blazem::norm(ma-*mb);
        ret+=norm_ab*norm_ab;
        mb++;
    }
    return sqrt(ret);
}

template <class T> double relative_fnorm(const Projector3<T>& a, const Projector3<T>& b)
{
    double ret=0.0;
    assert(a.dense.size()==b.dense.size());
    auto mb=b.dense.begin();
    for (auto ma:a.dense)
    {
        double norm_ab=blazem::norm(ma-*mb);
        double avg_norm_ab=(blazem::norm(ma)+blazem::norm(*mb))/2.0;
        if (avg_norm_ab>0.0) norm_ab/=avg_norm_ab;
        ret+=norm_ab*norm_ab;
        mb++;
    }
    return sqrt(ret);
}

template double          fnorm<double>(const Projector3<double>& a, const Projector3<double>& b);
template double relative_fnorm<double>(const Projector3<double>& a, const Projector3<double>& b);

} // namespace qchem
