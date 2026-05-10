// File: BasisSet/Atom/radial/BSpline/Imp/SplineGrouper.C
module;
#include <algorithm>
#include <iostream>
module qchem.BasisSet.Atom.BSpline.SplineGrouper;
import qchem.Types;
using std::cout;
using std::endl;

template <size_t K> size_t SplineGrouper<K>::Insert(const spline_t& sp,size_t l)
{
    // double rmin=sp.getSupport().front();
    size_t index=unique_sp.size();
    if (const auto &ie =unique_sp.find(sp);ie==unique_sp.end())
    {
        // new exponent
        unique_spv.push_back(sp);
        maxls.push_back(l);
        unique_sp[sp]=index;
    }
    else 
        index=ie->second;

    if (l>maxls[index]) maxls[index]=l;
    // cout << "SplineGrouper index,rmin,l,maxl=" << index << " " << rmin << " " << l << " " << maxls[index] << endl;
    return index;
}
 
#define INSTANCEk(k) template class SplineGrouper<k>;
#include "../Instance.hpp"

