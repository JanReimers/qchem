// File: BasisSetImp/TCommon.H

#include "Imp/BasisSet/TCommon.H"
#include "Imp/Containers/ptr_vector_io.h"
#include <Irrep_BS.H>

void BasisSetImp::Set(const LAParams& lap)
{
    for (auto b:itsBasisSets) b->Set(lap);
}

size_t BasisSetImp::GetNumFunctions() const
{
    index_t ret=0;
    for (auto bs:*this) 
        ret+=bs->GetNumFunctions();
    return ret;
}

void BasisSetImp::Insert(bs_t* bs)
{
    assert(bs);
    itsBasisSets.push_back(bs);
}


//
//  StreamableObject stuff.
//
std::ostream&  BasisSetImp::Write(std::ostream&  os  ) const
{
    return os << itsBasisSets;
}

std::istream&  BasisSetImp::Read (std::istream&  is  )
{
    return is;
}
