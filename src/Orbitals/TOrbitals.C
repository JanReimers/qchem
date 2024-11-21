// File: OrbitalGroupImplementation.C  general orbital group implementation.



#include "Imp/Orbitals/TOrbitals.H"
#include "Imp/Orbitals/TOrbital.H"
#include "Imp/ChargeDensity/IrrepCD.H"
#include "Imp/Misc/DFTDefines.H"
#include <Hamiltonian.H>
#include <QuantumNumber.H>
#include <LASolver.H>
#include "Imp/Containers/ptr_vector_io.h"
#include "oml/vector.h"
#include "oml/smatrix.h"
#include "oml/matrix.h"
#include <iostream>

//-----------------------------------------------------------------
//
//  Construction zone
//
template <class T> TOrbitalsImp<T>::TOrbitalsImp()
{};

template <class T> TOrbitalsImp<T>::
TOrbitalsImp(const TIrrepBasisSet<T>* bs)
    : itsBasisSet(bs)
    , itsLASolver(bs->CreateSolver())
{};


//-----------------------------------------------------------------
//
//  Orbitals stuff.
//
template <class T> index_t TOrbitalsImp<T>::GetNumOrbitals() const
{
    return itsOrbitals.size();
}

template <class T> double TOrbitalsImp<T>::GetEigenValueChange(const Orbitals& og) const
{
    // No UT coverage
    // TODO: OrbitalGroup should return a vector of energies.
    double del=0;
    auto b2=og.begin();
    for (auto b1:*this)
    {
        del+=Square(b1->GetEigenEnergy()-(*b2)->GetEigenEnergy());
        b2++;
    }
    return sqrt(del);
}

//
//  This is where the real SCF work gets done.
//
template <class T> void TOrbitalsImp<T>::UpdateOrbitals(const Hamiltonian& ham,const Spin& spin)
{
    assert(itsBasisSet);
    SMatrix<T> H=ham.BuildHamiltonian(itsBasisSet,spin);
//    std::cout << "UpdateOrbitals " << itsBasisSet->GetQuantumNumber() << " spin=" << spin << std::endl;
    assert(!isnan(H));
    auto [U,e]=itsLASolver->Solve(H);
    itsOrbitals.clear();
    index_t n=e.size();
    //std::cout << "   Eigen values=" << e << std::endl;
    for (index_t i=1; i<=n; i++)
        itsOrbitals.push_back(new TOrbitalImp<T>(itsBasisSet,U.GetColumn(i), e(i),spin, i));
}


template <class T> ChargeDensity* TOrbitalsImp<T>::GetChargeDensity(Spin s) const
{
    return new IrrepCD<T>(CalculateDensityMatrix(),itsBasisSet,s);
}


template <class T> typename TOrbitalsImp<T>::SMat TOrbitalsImp<T>::
CalculateDensityMatrix() const
{
    SMat d(itsBasisSet->GetNumFunctions());
    Fill(d,T(0.0));
    for (auto b=this->beginT();b!=this->end();b++) b->AddDensityMatrix(d);
    return d;
}

//-----------------------------------------------------------------
//
//  VectorFunction stuff.
//
template <class T> typename TOrbitalsImp<T>::Vec TOrbitalsImp<T>::
operator()(const RVec3& r) const
{
    Vec ret(GetNumOrbitals());
    typename Vec::iterator i(ret.begin());
    // No UT coverage
    for (auto b=this->beginT();b!=this->end();i++,b++) *i=(**b)(r);
    return ret;
}

template <class T> typename TOrbitalsImp<T>::Vec3Vec TOrbitalsImp<T>::
Gradient(const RVec3& r) const
{
    // No UT coverage
    Vec3Vec ret(GetNumOrbitals());
    typename Vec3Vec::iterator i(ret.begin());
    for (auto b=this->beginT();b!=this->end();i++,b++) *i=b->Gradient(r);
    return ret;
}

//-----------------------------------------------------------------
//
//  Streamable stuff.
//
template <class T> std::ostream& TOrbitalsImp<T>::Write(std::ostream& os) const
{
     if (!StreamableObject::Pretty())
    {
        os  << itsOrbitals;
        if (!StreamableObject::Binary()) os << std::endl;
    }
    else
    {
        os << "        Orbital group with " << GetNumOrbitals() << " " << itsBasisSet->GetQuantumNumber() << "orbitals:" << std::endl;
        os << "            Occupation      Energy      Eigenvector" << std::endl;
        os << itsOrbitals;
    }
    if (!StreamableObject::Pretty())
    {
        os  << itsBasisSet;
        if (StreamableObject::Ascii()) os << std::endl;
    }
    return os;
}

template <class T> std::istream& TOrbitalsImp<T>::Read(std::istream& is)
{
    is >> itsOrbitals;
    if (!StreamableObject::Binary()) is >> std::ws;
    return is;
}


template class TOrbitalsImp<double>;
//template class TOrbitalGroupImplementation<std::complex<double> >;
