// File: BasisSet.C  Basis set factory.

#include "BasisSet.H"
#include <complex>

template class TBasisFunction<double>;
template class TBasisFunction<std::complex<double> >;
template class TIrrepBasisSet<double>;
template class TIrrepBasisSet<std::complex<double> >;

#include "DFTDataBase/HeapDB/HeapDB.H"
#include "Imp/Containers/ptr_vector_io.h"

BasisGroup::BasisGroup()
: itsDB(new HeapDB<double>())
, itsBasisSets()
{
}

BasisGroup::BasisGroup(AnalyticIE<double>* ie)
: itsDB(new HeapDB<double>(ie,this))
, itsBasisSets()
{
}

BasisGroup::~BasisGroup() 
{
    delete itsDB;
};

size_t BasisGroup::GetNumFunctions() const
{
    index_t ret=0;
    for (auto bs:*this) ret+=bs->GetNumFunctions();
    return ret;
}

size_t BasisGroup::GetNumBasisSets() const
{
    return itsBasisSets.size();
}

void BasisGroup::Insert(IrrepBasisSet* bs)
{
    assert(bs);
    bs->SetStartIndex(GetNumFunctions()+1);
    itsBasisSets.push_back(bs);
}

BasisGroup::iecv_t BasisGroup::Flatten() const
{
    iecv_t iecs;
    for (auto bs:*this)
    {
        //const TIrrepBasisSet<double>* tbs=dynamic_cast<const TIrrepBasisSet<double>*>(bs); //TODO we need a way to avoid all the TBasisSet casts
        iecs.push_back(bs);
    }
    return iecs;
}

//
//  StreamableObject stuff.
//
std::ostream&  BasisGroup::Write(std::ostream&  os  ) const
{
    return os << itsBasisSets;
}

std::istream&  BasisGroup::Read (std::istream&  is  )
{
    return is;
}
