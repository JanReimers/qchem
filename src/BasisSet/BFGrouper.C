// File: BFGrouper.C  Group Slater or Gaussian basis functions by unique exponents.

#include "Imp/BasisSet/BFGrouper.H"
#include <cassert>

void BFGrouper::Append(double e, size_t l, size_t flat_index)
{
    size_t index=unique_es.size();
    if (const auto &ie =unique_es.find(e);ie==unique_es.end())
    {
        unique_esv.push_back(e);
        unique_es[e]=index;
    }
    else 
        index=ie->second;

    es_indices.push_back(index);  
    
    if (const auto &il =L_indices.find(l);il==L_indices.end())
        L_indices[l]=std::vector<size_t>();

    L_indices[l].push_back(flat_index);
}

const std::vector<size_t>& BFGrouper::indices(size_t l) const
{
    auto i=L_indices.find(l);
    assert(i!=L_indices.end());
    return i->second;
}

