// File: BasisSetImp/TCommon.H

#include "Imp/BasisSet/BS_Common.H"
#include "Common/stl_io.h"
#include <Irrep_BS.H>

void BS_Common::Insert(bs_t* bs)
{
    assert(bs);
    itsBasisSets.push_back(std::unique_ptr<bs_t>(bs));
}

void BS_Common::Set(const LAParams& lap)
{
    for (auto& b:itsBasisSets) b->Set(lap);
}

size_t BS_Common::GetNumFunctions() const
{
    index_t ret=0;
    for (auto& bs:*this) 
        ret+=bs->GetNumFunctions();
    return ret;
}

BasisSet::symv_t BS_Common::GetSymmetries  () const
{
    symv_t symv;
    for (auto& b:itsBasisSets) symv.push_back(b->GetSymmetry());
    return symv;
}




//
//  StreamableObject stuff.
//
std::ostream&  BS_Common::Write(std::ostream&  os  ) const
{
    return os << itsBasisSets;
}

