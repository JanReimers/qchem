// File: NumericalIE.C  Here is where all the numerical integrals get calculated.



#include "BasisSetImplementation/NumericalIEImp.H"
#include "LASolver/LASolver.H"
#include "SCFIterator/IterationParams.H"

#include "BasisSet.H"
#include "Cluster.H"
#include "Mesh/Mesh.H"
#include "oml/smatrix.h"
#include "oml/imp/binio.h"
#include <iostream>
#include <cassert>
#include <cstdlib>

#define TYPE_STRING "TBasisSet<double>"
#define TYPE TBasisSet<double>
#include "Misc/Persistent/IDRef.Ci"


template <class T> typename NumericalIE<T>::RSMat NumericalIE<T>::
    MakeInverse(const RSMat& S) 
{
    LinearAlgebraParams lap={qchem::Lapack,qchem::SVD,1e-6,1e-12};
    LASolver<double>* las=LASolver<double>::Factory(lap);
    return las->Inverse(S);
}

template typename NumericalIE<double>::RSMat NumericalIE<double>::MakeInverse(const RSMat& S) ;
//-----------------------------------------------------------------
//
//  Construction zone.
//

template <class T> NumericalIEImp<T>::NumericalIEImp()
    : itsMesh                        (0)
    , itsIntegrator                  (0)
{}; //No UT coverage.

template <class T> NumericalIEImp<T>::NumericalIEImp(const NumericalIEImp& ie)
    : itsMesh         (ie.itsMesh->Clone())
    , itsIntegrator(new MeshIntegrator<T>(itsMesh))
{};//No UT coverage.

template <class T> NumericalIEImp<T>::NumericalIEImp(Mesh* theMesh)
    : itsMesh               (theMesh    )
    , itsIntegrator(new MeshIntegrator<T>(itsMesh))
{};

template <class T> NumericalIEImp<T>::~NumericalIEImp()
{
    delete itsIntegrator; //Integrator does not own the mesh.
    delete itsMesh;
}

template <class T> void NumericalIEImp<T>::Insert(const TBasisSet<T>* theSet)
{
    itsBasisSet=IDRef<const TBasisSet<T> >(theSet);
//    itsNormalizations.SetLimits(VecLimits(itsN));
//    InitNormalizations();    
    itsNormalizations = GetNumericalNormalization(*itsBasisSet);
    assert(!isnan(itsNormalizations));
    CheckInitialized();
}

template <class T> const typename NumericalIEImp<T>::RVec NumericalIEImp<T>::GetNumericalNormalization(const TBasisSet<T>& bs) const
{
    CheckInitialized();
    RVec ret=itsIntegrator->Normalize(bs);
    return ret;
}

template <class T> void NumericalIEImp<T>::CheckInitialized() const
{
}

//-----------------------------------------------------------------
//
//  Streamable Object stuff
//
template <class T> std::ostream& NumericalIEImp<T>::Write(std::ostream& os) const
{
    os << itsBasisSet << itsNormalizations;
    return os << itsMesh;
}

template <class T> std::istream& NumericalIEImp<T>::Read (std::istream& is)
{
    is >> itsBasisSet >> itsNormalizations;
    delete itsMesh;
    itsMesh=Mesh::Factory(is);
    is >> itsMesh;
    delete itsIntegrator;
    itsIntegrator=new MeshIntegrator<T>(itsMesh);
    return is;
}

template <class T> NumericalIE<T>* NumericalIEImp<T>::Clone() const
{
    return new NumericalIEImp(*this);
}

template <class T> void NumericalIEImp<T>::
Normalize(const RVec& n1, Mat& m, const RVec& n2) const
{
    typename  Mat::Subscriptor      s(m);
    for (unsigned int i=1; i<=n1.size(); i++)
        for (unsigned int j=1; j<=n2.size(); j++)
            s(i,j)*=n1(i)*n2(j);
}

template <class T> typename NumericalIEImp<T>::RVec NumericalIEImp<T>::MakeNormalization() const
{
    // No UT coverage
    CheckInitialized();
    RVec ret=itsNormalizations;
    int i=1;
    for (auto b:*itsBasisSet) ret(i++)=b->GetNormalization();
    return ret;
}

template <class T> typename NumericalIEImp<T>::RVec NumericalIEImp<T>::MakeCharge() const
{
    //No UT coverage.
    CheckInitialized();
    RVec ret=itsIntegrator->Integrate(*itsBasisSet);
    ret=DirectMultiply(ret,itsNormalizations);
    return ret;
}

//--------------------------------------------------------------------------------------------
//
//  Overlap type integrals.
//
template <class T> typename NumericalIEImp<T>::SMat NumericalIEImp<T>::MakeOverlap() const
{
    //No UT coverage.
    CheckInitialized();
    std::cout << "Doing numerical overlap ...";
    std::cout.flush();
    SMat ret=itsIntegrator->Overlap(*itsBasisSet);
    Normalize(ret);
    std::cout << " done" << std::endl;
    return ret;
}

template <class T> typename NumericalIEImp<T>::Mat NumericalIEImp<T>::MakeOverlap(const TBasisSet<T>& theOtherBasisSet) const
{
    //No UT coverage.
    CheckInitialized();
//
//  Can't assume other basis set is numerical, so we must explicitly evaluate
//  the normalization constants over this mesh.
//
    RVec otherNormalizations = GetNumericalNormalization(theOtherBasisSet);
    Mat ret=itsIntegrator->Overlap(*itsBasisSet,theOtherBasisSet);
    Normalize(itsNormalizations,ret,otherNormalizations);
    return ret;
}



template <class T> typename NumericalIEImp<T>::Vec NumericalIEImp<T>::MakeOverlap(const ScalarFunction<double>& f) const
{
    CheckInitialized();
    Vec ret=itsIntegrator->Overlap(f,*itsBasisSet);
    assert(!isnan(ret));
    Normalize(ret);
    assert(!isnan(ret));
    return ret;
}

//--------------------------------------------------------------------------------------------
//
//  Repulsion type integrals are not implementated, they are 6 dimensional!
//
template <class T> typename NumericalIEImp<T>::SMat NumericalIEImp<T>::MakeRepulsion() const
{
    //No UT coverage.
    std::cerr << "Error: NumericalIE<T>::MakeRepulsion Do not do repulsion integrals numerically" << std::endl;
    assert(false);
    CheckInitialized();
    SMat ret=itsIntegrator->Repulsion(*itsBasisSet);
    assert(!isnan(ret));
    Normalize(ret);
    assert(!isnan(ret));
    return ret;
}

template <class T> typename NumericalIEImp<T>::Mat NumericalIEImp<T>::MakeRepulsion(const TBasisSet<T>& theOtherBasisSet) const
{
    //No UT coverage.
    std::cerr << "Error: NumericalIE<T>::MakeRepulsion Do not do repulsion integrals numerically" << std::endl;
    assert(false);
    CheckInitialized();
//
//  Can't assume other basis set is numerical, so we must explicitly evaluate
//  the normalization constants over this mesh.
//
    RVec otherNormalizations = GetNumericalNormalization(theOtherBasisSet);
    Mat ret=itsIntegrator->Repulsion(*itsBasisSet,theOtherBasisSet);
    Normalize(itsNormalizations,ret,otherNormalizations);
    return ret;
}

template <class T> typename NumericalIEImp<T>::Vec NumericalIEImp<T>::MakeRepulsion(const ScalarFunction<double>& f) const
{
    //No UT coverage.
    CheckInitialized();
    Vec ret=itsIntegrator->Repulsion(f,*itsBasisSet);
    Normalize(ret);
    return ret;
}


#ifdef USE_FOR_DEBUGGING_ANALYTIC

template <class T> void NumericalIEImp<T>::MakeOverlap3C(ERI3& mlist, const TBasisSet<T>& bs) const
{
    // No UT coverage
    for (auto b=bs.beginT(); b!=bs.end(); b++) mlist.Add(MakeOverlap(**b));
}

template <class T> void NumericalIEImp<T>::MakeRepulsion3C(ERI3& ret,const TBasisSet<T>& bs) const
{
    //No UT coverage.
    std::cerr << "Error: NumericalIE<T>::MakeRepulsion Do not do repulsion integrals numerically" << std::endl;
    assert(false);
    // No UT coverage
    for (auto b=bs.beginT(); b!=bs.end(); b++) ret.Add(MakeRepulsion(**b));
}

template <class T> void NumericalIEImp<T>::MakeRepulsion4C(ERIList& eris,ERIList&,const BasisGroup*) const
{
    //No UT coverage.
    std::cerr << "NumericalIE<T>::MakeRepulsion4C 4 center numerical integrals are not supported" << std::endl;
    exit(-1);
}
template <class T> typename NumericalIEImp<T>::SMat NumericalIEImp<T>::MakeKinetic() const
{
    //No UT coverage.
    CheckInitialized();
    SMat ret=itsIntegrator->Grad(*itsBasisSet);

    Normalize(ret);
    ret  *=  0.5;
    return ret;
}

template <class T> typename NumericalIEImp<T>::SMat NumericalIEImp<T>::MakeNuclear(const Cluster& theCluster) const
{
    //No UT coverage.
    CheckInitialized();
    SMat ret=itsIntegrator->Nuclear(*itsBasisSet);

    Normalize(ret);
    ret  *=  -(double)(theCluster.GetNuclearCharge());
    return ret;
}


#endif


//----------------------------------------------------------
//
//  Private stuff.
//
template <class T> typename NumericalIEImp<T>::SMat NumericalIEImp<T>::MakeOverlap(const TBasisFunction<T>& bf) const
{
    //No UT coverage.
    CheckInitialized();
    SMat ret=itsIntegrator->Overlap3C(*itsBasisSet,bf);
    Normalize(ret);
    return ret;
}

template <class T> typename NumericalIEImp<T>::SMat NumericalIEImp<T>::MakeRepulsion(const TBasisFunction<T>& bf) const
{
    //No UT coverage.
    CheckInitialized();
    SMat ret = itsIntegrator->Repulsion3C(*itsBasisSet,bf);
    Normalize(ret);
    return ret;
}


template class NumericalIEImp<double>;
template class NumericalIEImp<std::complex<double> >;


