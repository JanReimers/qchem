// File: OrbitalGroupImplementation.C  general orbital group implementation.



#include "Imp/Orbitals/TOrbitals.H"
#include "Imp/Orbitals/TOrbital.H"
#include "Imp/ChargeDensity/IrrepCD.H"
#include "Imp/Misc/DFTDefines.H"
#include "Imp/SCFAccelerator.H"
#include <Irrep_BS.H>
#include <Hamiltonian.H>
#include <Symmetry.H>
#include <Orbital_QNs.H>
#include <LASolver.H>
#include "Imp/Containers/stl_io.h"
#include "oml/vector.h"
#include "oml/smatrix.h"
#include "oml/matrix.h"
#include <iostream>

//-----------------------------------------------------------------
//
//  Construction zone
//
template <class T> TOrbitalsImp<T>::
TOrbitalsImp(const TOrbital_IBS<T>* bs, Spin ms,SCFIrrepAccelerator* acc)
    : itsBasisSet(bs)
    , itsLASolver(bs->CreateSolver())
    , itsQNs(ms,&bs->GetSymmetry())
    , itsAccelerator(acc)
{
    assert(itsBasisSet->GetNumFunctions()>0);
    assert(itsAccelerator);
    itsAccelerator->Init(itsLASolver);
};

template <class T> TOrbitalsImp<T>::~TOrbitalsImp()
{
    delete itsLASolver;
    delete itsAccelerator;
}

//-----------------------------------------------------------------
//
//  Orbitals stuff.
//
template <class T> index_t TOrbitalsImp<T>::GetNumOrbitals() const
{
    return itsOrbitals.size();
}
template <class T> index_t TOrbitalsImp<T>::GetNumOccOrbitals() const
{
    index_t n=0;
    for (auto& o:*this)
        if (o->IsOccupied()) 
            n++;

    return n;
}
template <class T> double TOrbitalsImp<T>::GetEigenValueChange(const Orbitals& og) const
{
    // No UT coverage
    // TODO: OrbitalGroup should return a vector of energies.
    double del=0;
    auto b2=og.Iterate<Orbital>().begin();
    for (auto b1:Iterate<Orbital>())
    {
        const Orbital* o2=*b2;
        del+=Square(b1->GetEigenEnergy()-o2->GetEigenEnergy());
        ++b2;
    }
    return sqrt(del);
}

//
//  This is where the real SCF work gets done.
//
template <class T> void TOrbitalsImp<T>::UpdateOrbitals(Hamiltonian& ham,const DM_CD* cd)
{
    assert(itsBasisSet);
    SMatrix<T> H=ham.GetMatrix(itsBasisSet,itsQNs.ms,cd);
    //std::cout << "UpdateOrbitals " << itsBasisSet->GetSymmetry() << " spin=" << spin << std::endl;
    //std::cout << "H=" << H << std::endl;
    assert(!isnan(H));
    auto [U,e]=itsLASolver->Solve(H);

    itsOrbitals.clear();
    index_t n=e.size();
    //std::cout << "   Eigen values=" << e << std::endl;

    //
    //  Strip out all the positron orbitals.
    //
    static const double e_positron=-c_light*c_light; //Max positron energy = -2mc^2.
    size_t index=1;
    for (index_t i=1; i<=n; i++)
    {
 //               std::cout << "o=" << o->GetEigenEnergy() << std::endl;
        if (e(i)<=e_positron) continue; //Strip out all the positron orbitals.
        Orbital* o=new TOrbitalImp<T>(itsBasisSet,U.GetColumn(i), e(i),Orbital_QNs(index++,itsQNs));
        itsOrbitals.push_back(std::unique_ptr<Orbital>(o));

    }
}


template <class T> DM_CD* TOrbitalsImp<T>::GetChargeDensity() const
{
    return new IrrepCD<T>(CalculateDensityMatrix(),itsBasisSet,GetQNs());
}


template <class T> typename TOrbitalsImp<T>::SMat TOrbitalsImp<T>::
CalculateDensityMatrix() const
{
    SMat d(itsBasisSet->GetNumFunctions());
    Fill(d,T(0.0));
    for (auto b:Iterate<TOrbital<double>>()) b->AddDensityMatrix(d);
    return d;
}

template <class T>  Irrep_QNs TOrbitalsImp<T>::GetQNs() const
{
    return itsQNs;
}

template <class T>  typename TOrbitalsImp<T>::RVec TOrbitalsImp<T>::Get_BS_Diagonal() const
{
    assert(itsLASolver);
    return itsLASolver->Get_BS_Diagonal();
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
    for (auto b:Iterate<TOrbital<double>>()) 
    {
        *i=(*b)(r);
        i++;
    }
    return ret;
}

template <class T> typename TOrbitalsImp<T>::Vec3Vec TOrbitalsImp<T>::
Gradient(const RVec3& r) const
{
    // No UT coverage
    Vec3Vec ret(GetNumOrbitals());
    typename Vec3Vec::iterator i(ret.begin());
    for (auto b:Iterate<TOrbital<double>>()) 
    {
        *i=b->Gradient(r);
        i++;
    } 
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
        os << "        Orbital group with " << GetNumOrbitals() << " " << itsBasisSet->GetSymmetry() << "orbitals:" << std::endl;
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


template class TOrbitalsImp<double>;
