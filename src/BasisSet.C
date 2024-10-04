// File: BasisSet.C  Basis set factory.

#include "BasisSet.H"
#include <complex>

template class TBasisFunction<double>;
template class TBasisFunction<std::complex<double> >;
template class TBasisSet<double>;
template class TBasisSet<std::complex<double> >;

#include "IntegralDataBase.H"
#include "IntegralEngine.H"
#include "Misc/ptr_vector1_io.h"

BasisGroup::BasisGroup()
: itsBasisSets()
{
}

BasisGroup::~BasisGroup() {};

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

BasisGroup::iev_t BasisGroup::Flatten() const
{
    iev_t ies;
    for (auto bs:*this)
    {
        const TBasisSet<double>* tbs=dynamic_cast<const TBasisSet<double>*>(bs); //TODO we need a way to avoid all the TBasisSet casts
        ies.push_back(tbs->GetIntegralEngine1());
    }
    return ies;
}

void BasisGroup::Insert(const ERIList& C, const ERIList& X) const
{
    for (auto bs:*this)
    {
        const TBasisSet<double>* tbs=dynamic_cast<const TBasisSet<double>*>(bs); //TODO we need a way to avoid all the TBasisSet casts
        tbs->GetDataBase()->Insert(C,X);
    }
}
void BasisGroup::Insert(const ERIList1& J, const ERIList1& K) const
{
    for (auto bs:*this)
    {
        const TBasisSet<double>* tbs=dynamic_cast<const TBasisSet<double>*>(bs); //TODO we need a way to avoid all the TBasisSet casts
        tbs->GetDataBase()->Insert(J,K);
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
