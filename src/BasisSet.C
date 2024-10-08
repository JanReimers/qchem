// File: BasisSet.C  Basis set factory.

#include "BasisSet.H"
#include <complex>

template class TBasisFunction<double>;
template class TBasisFunction<std::complex<double> >;
template class TIrrepBasisSet<double>;
template class TIrrepBasisSet<std::complex<double> >;

#include "IntegralDataBase.H"
#include "NumericalIE.H"
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

void BasisGroup::Insert(IrrepBasisSet* bs)
{
    assert(bs);
    int N=GetNumFunctions();
    bs->SetStartIndex(N+1);
    bs->Insert(this);
    itsBasisSets.push_back(bs);
    TIrrepBasisSet<double>* tbs=dynamic_cast<TIrrepBasisSet<double>*>(bs);
    assert(tbs);
    tbs->GetDataBase()->Insert(this);
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

void BasisGroup::Insert(const ERI4& J, const ERI4& K) const
{
    for (auto bs:*this)
    {
        const TIrrepBasisSet<double>* tbs=dynamic_cast<const TIrrepBasisSet<double>*>(bs); //TODO we need a way to avoid all the TBasisSet casts
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
