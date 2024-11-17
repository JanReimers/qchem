// File: BFGrouper.C  Group Slater or Gaussian basis functions by unique exponents.

#include "Imp/BasisSet/BFGrouper.H"
#include <cassert>

void BFGrouper::Append(double e, size_t l)
{
    size_t index=unique_es.size();
    if (const auto &ie =unique_es.find(e);ie==unique_es.end())
    {
        unique_esv.push_back(e);
        maxls.push_back(l);
        unique_es[e]=index;
    }
    else 
        index=ie->second;

    if (l>maxls[index]) maxls[index]=l;
    es_indices.push_back(index);  
}

//  indices should already be zero based.
size_t BFGrouper::LMax(size_t ia, size_t ib, size_t ic, size_t id) const
{
    size_t lmax_ab=std::max(maxls[ia],maxls[ib]);
    size_t lmax_cd=std::max(maxls[ic],maxls[id]);
    return std::max(lmax_ab,lmax_cd);
}


