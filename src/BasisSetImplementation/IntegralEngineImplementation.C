// File: IntegralEngineImplementation.C  Implement common attributes for Integral Engines.



#include "BasisSetImplementation/IntegralEngineImplementation.H"
#include "BasisSet/TBasisSet.H"
#include <iostream>
#include <cassert>

#define TYPE_STRING "TBasisSet"
#define TYPE TBasisSet<double>
#include "Misc/Persistent/IDRef.Ci"

template class IDRef<const TBasisSet<double> >;
template class IDRef<const TBasisSet<std::complex<double> > >;

//-----------------------------------------------------------------
//
//  Construction zone.
//

template <class T> IntegralEngineImplementation<T>::IntegralEngineImplementation()
    : itsBasisSet      ( )
    , itsN             (0)
    , itsNormalizations( )
{};

template <class T> IntegralEngineImplementation<T>::IntegralEngineImplementation(const IntegralEngineImplementation&)
    : itsBasisSet      ( )
    , itsN             (0)
    , itsNormalizations( )
{};

template <class T> void IntegralEngineImplementation<T>::Insert(const TBasisSet<T>* theBasisSet)
{
    itsBasisSet=IDRef<const TBasisSet<T> >(theBasisSet);
    itsN=theBasisSet->GetNumFunctions();
    itsNormalizations.SetLimits(VecLimits(itsN));
    InitNormalizations();
    CheckBasisSet();
}

//----------------------------------------------------------------------------------------
//
//  Check that initialization was doen properly.
//
template <class T> void IntegralEngineImplementation<T>::CheckBasisSet() const
{
    assert(&*itsBasisSet);
    assert(itsN==itsBasisSet->GetNumFunctions());
    assert(itsN==itsNormalizations.size());
    assert(!isnan(itsNormalizations));
}

//----------------------------------------------------------------------------------------
//
//  Streamable stuff.
//
template <class T> std::ostream& IntegralEngineImplementation<T>::Write(std::ostream& os) const
{
    assert(!isnan(itsNormalizations));
    return os << itsBasisSet << itsNormalizations;
}

template <class T> std::istream& IntegralEngineImplementation<T>::Read (std::istream& is)
{
    is >> itsBasisSet;
    is >> itsNormalizations;
    assert(!isnan(itsNormalizations));
    return is;
}

//----------------------------------------------------------------------------------------
//
//  Initialize all the local info.
//
template <class T> void IntegralEngineImplementation<T>::InitNormalizations()
{
    assert(itsN==itsNormalizations.size());
    int i=1;
    for (auto b:*itsBasisSet) itsNormalizations(i++)=b->GetNormalization();
    assert(!isnan(itsNormalizations));
}

//----------------------------------------------------------------------------------------
//
//  Special integrals
//
template <class T> typename IntegralEngineImplementation<T>::RVec IntegralEngineImplementation<T>::
MakeNormalization() const
{
    CheckBasisSet();
    return itsNormalizations;
}

template <class T> typename IntegralEngineImplementation<T>::RVec IntegralEngineImplementation<T>::
MakeCharge() const
{
    CheckBasisSet();
    RVec ret(itsN);
    int i=1;
    for (auto b:*itsBasisSet) ret(i++)=b->GetCharge();
    return DirectMultiply(ret,itsNormalizations);
}

template <class T> void IntegralEngineImplementation<T>::
Normalize(const RVec& n1, Mat& m, const RVec& n2) const
{
    typename  Mat::Subscriptor      s(m);
    for (index_t i=1; i<=n1.size(); i++)
        for (index_t j=1; j<=n2.size(); j++)
            s(i,j)*=n1(i)*n2(j);
}

template <class T> std::vector<typename IntegralEngineImplementation<T>::SMat> IntegralEngineImplementation<T>::MakeMatrixList(int n) const
{
    std::vector<SMat> ret;
    for (index_t i=0; i<n; i++)
    {
        ret.push_back(SMat(itsN,itsN));
        Fill(ret.back(),T(0.0));
    }
    return ret;
}

template class IntegralEngineImplementation<double>;
template class IntegralEngineImplementation<std::complex<double> >;
