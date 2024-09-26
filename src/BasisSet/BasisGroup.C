// File BasisGroup.C  Store a set of basis functions and manage the ERI super matrix.

#include "BasisSet/BasisGroup.H"
#include "BasisSet/BasisGroupBrowser.H"
#include "BasisSet/TBasisSet.H"
#include "BasisSet/IntegralDataBase.H"
#include "BasisSet/IntegralEngine.H"
#include "Misc/ptr_vector_io.h"

BasisGroup::BasisGroup()
: itsBasisSets()
{
}

BasisGroup::~BasisGroup() {};

index_t BasisGroup::GetNumFunctions() const
{
    index_t ret=0;
    for (BasisGroupBrowser b(*this);b;b++) ret+=b->GetNumFunctions();
    return ret;
}

index_t BasisGroup::GetNumBasisSets() const
{
    return itsBasisSets.size();
}

void BasisGroup::Insert(BasisSet* bs)
{
    assert(bs);
    int N=GetNumFunctions();
    bs->SetStartIndex(N+1);
    itsBasisSets.push_back(bs);
    TBasisSet<double>* tbs=dynamic_cast<TBasisSet<double>*>(bs);
    assert(tbs);
    tbs->GetDataBase()->Insert(this);
}

void BasisGroup::Insert(const ERIList& C, const ERIList& X) const
{
    for (BasisGroupBrowser bs(*this);bs;bs++)
    {
        const TBasisSet<double>* tbs=dynamic_cast<const TBasisSet<double>*>(&bs); //TODO we need a way to avoid all the TBasisSet casts
        tbs->GetDataBase()->Insert(C,X);
    }
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
