// File: BasisSet/Atom/radial/BSpline/Imp/SplineGrouper.C
module;
#include <algorithm>
module qchem.BasisSet.Atom.Internal.SplineGrouper;
import qchem.Types;

template <size_t K> size_t SplineGrouper<K>::Insert(const spline_t& sp,size_t l)
{
    double rmin=sp.getSupport().front();
    size_t index=unique_es.size();
    if (const auto &ie =unique_es.find(rmin);ie==unique_es.end())
    {
        // new exponent
        unique_spv.push_back(sp);
        maxls.push_back(l);
        unique_es[rmin]=index;
    }
    else 
        index=ie->second;

    if (l>maxls[index]) maxls[index]=l;
    return index;
}
 
#define INSTANCEk(k) template class SplineGrouper<k>;
#include "../Instance.hpp"

