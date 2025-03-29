// File: BasisSetImp/TCommon.H

#include "Imp/BasisSet/BS_Common.H"
#include "Imp/Containers/ptr_vector_io.h"
#include <Irrep_BS.H>

void BS_Common::Set(const LAParams& lap)
{
    for (auto b:itsBasisSets) b->Set(lap);
}

size_t BS_Common::GetNumFunctions() const
{
    index_t ret=0;
    for (auto bs:*this) 
        ret+=bs->GetNumFunctions();
    return ret;
}

void BS_Common::Insert(bs_t* bs)
{
    assert(bs);
    itsBasisSets.push_back(bs);
}


//
//  StreamableObject stuff.
//
std::ostream&  BS_Common::Write(std::ostream&  os  ) const
{
    return os << itsBasisSets;
}

std::istream&  BS_Common::Read (std::istream&  is  )
{
    return is;
}
