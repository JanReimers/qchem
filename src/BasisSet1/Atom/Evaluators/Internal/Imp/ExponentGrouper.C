// File: BasisSet/Atom/radial/Imp/ExponentGrouper.C
module;
#include <algorithm>
module qchem.BasisSet.Atom.Internal.ExponentGrouper;
import qchem.Types;

size_t ExponentGrouper::Insert(double e,size_t l)
{
    size_t index=unique_es.size();
    if (const auto &ie =unique_es.find(e);ie==unique_es.end())
    {
        // new exponent
        unique_esv.push_back(e);
        maxls.push_back(l);
        unique_es[e]=index;
    }
    else 
        index=ie->second;

    if (l>maxls[index]) maxls[index]=l;
    return index;
}
 


//  indices should already be zero based.
size_t Grouper::LMax(size_t ia, size_t ib, size_t ic, size_t id) const
{
    size_t lmax_ab=std::max(maxls[ia],maxls[ib]);
    size_t lmax_cd=std::max(maxls[ic],maxls[id]);
    return std::max(lmax_ab,lmax_cd);
}


