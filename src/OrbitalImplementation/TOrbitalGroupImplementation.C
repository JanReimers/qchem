// File: OrbitalGroupImplementation.C  general orbital group implementation.



#include "OrbitalImplementation/TOrbitalGroupImplementation.H"
#include "OrbitalImplementation/TOrbitalImplementation.H"
#include "Orbital/TOrbitalGroupBrowser.H"
#include "ChargeDensityImplementation/ExactIrrepCD/ExactIrrepCD.H"
#include "oml/vector.h"
#include "oml/vector.h"
#include "oml/smatrix.h"
#include "oml/matrix.h"
#include <iostream>

//-----------------------------------------------------------------
//
//  Construction zone
//
template <class T> TOrbitalGroupImplementation<T>::TOrbitalGroupImplementation()
{};

template <class T> TOrbitalGroupImplementation<T>::
TOrbitalGroupImplementation(const rc_ptr<const BasisSet>& bs,
                            const Mat & evec,
                            const RVec& eval,
                            const Spin& S)
    : OrbitalGroupImplementation(bs)
{
    index_t n=eval.size();
    for (index_t i=1; i<=n; i++)
        itsOrbitals.push_back(new TOrbitalImplementation<T>(itsBasisSet,evec.GetColumn(i), eval(i),S));
};

//-----------------------------------------------------------------
//
//  VectorFunction stuff.
//
template <class T> typename TOrbitalGroupImplementation<T>::Vec TOrbitalGroupImplementation<T>::
operator()(const RVec3& r) const
{
    Vec ret(GetNumOrbitals());
    typename Vec::iterator i(ret.begin());
    TOrbitalGroupBrowser<T>  b(*this);
    for(; i!=ret.end()&&b; i++,b++)*i=(*b)(r);
    return ret;
}

template <class T> typename TOrbitalGroupImplementation<T>::Vec3Vec TOrbitalGroupImplementation<T>::
Gradient(const RVec3& r) const
{
    Vec3Vec ret(GetNumOrbitals());
    typename Vec3Vec::iterator i(ret.begin());
    TOrbitalGroupBrowser<T>  b(*this);
    for(; i!=ret.end()&&b; i++,b++) *i=b->Gradient(r);
    return ret;
}

//-----------------------------------------------------------------
//
//  Orbital stuff.
//
template <class T> ChargeDensity* TOrbitalGroupImplementation<T>::GetChargeDensity(Spin s) const
{
    return new ExactIrrepCD<T>(CalculateDensityMatrix(),itsRCBasisSet,s);
}

#include "BasisSet/QuantumNumber.H"
#include <iomanip>
template <class T> typename TOrbitalGroupImplementation<T>::SMat TOrbitalGroupImplementation<T>::
CalculateDensityMatrix() const
{
    index_t n=itsBasisSet->GetNumFunctions();
    SMat d(n,n);
    Fill(d,T(0.0));
    TOrbitalGroupBrowser<T> b(*this);
    std::cout.precision(4);
    for(; b; b++) b->AddDensityMatrix(d);
//	std::cout << "OrbitalGroup L=" << itsRCBasisSet->GetQuantumNumber() << "DensityMatrix=" << std::setw(7) << d << std::endl;
    return d;
}

//-----------------------------------------------------------------
//
//  Streamable stuff.
//
template <class T> std::ostream& TOrbitalGroupImplementation<T>::Write(std::ostream& os) const
{
    OrbitalGroupImplementation::Write(os);
    return os;
}

template <class T> std::istream& TOrbitalGroupImplementation<T>::Read(std::istream& is)
{
    OrbitalGroupImplementation::Read(is);
    return is;
}


template class TOrbitalGroupImplementation<double>;
template class TOrbitalGroupImplementation<std::complex<double> >;
