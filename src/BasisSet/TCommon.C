// File: BasisSetImp/TCommon.H

#include "Imp/BasisSet/TCommon.H"

#include "Imp/DataBase/HeapDB.H"
#include "Imp/Containers/ptr_vector_io.h"

BasisSetImp::BasisSetImp()
: itsDB(new HeapDB<double>())
, itsBasisSets()
{
}

BasisSetImp::BasisSetImp(AnalyticIE<double>* ie)
: itsIE(ie)
, itsDB(new HeapDB<double>(ie))
, itsBasisSets()
{
}

BasisSetImp::~BasisSetImp() 
{
    delete itsDB;
};

size_t BasisSetImp::GetNumFunctions() const
{
    index_t ret=0;
    for (auto bs:*this) ret+=bs->GetNumFunctions();
    return ret;
}

size_t BasisSetImp::GetNumIrreps() const
{
    return itsBasisSets.size();
}

void BasisSetImp::Insert(IrrepBasisSet* bs)
{
    assert(bs);
    bs->SetStartIndex(GetNumFunctions()+1);
    itsBasisSets.push_back(bs);
    bs->SetIndex(itsBasisSets.size());
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
