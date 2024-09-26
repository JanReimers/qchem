// File: NumericalIE.C  Here is where all the numerical integrals get calculated.



#include "BasisSetImplementation/NumericalIE.H"
#include "BasisSet/TBasisSet.H"
#include "BasisSet/TBasisFunction.H"
#include "BasisSet/TBasisSetBrowser.H"
#include "Cluster/Cluster.H"
#include "Mesh/Mesh.H"
#include "Misc/MatrixList.H"
#include "oml/imp/binio.h"
#include <iostream>
#include <cassert>
#include <cstdlib>

//-----------------------------------------------------------------
//
//  Construction zone.
//

template <class T> NumericalIE<T>::NumericalIE()
    : IntegralEngineImplementation<T>( )
    , itsMesh                        (0)
    , itsIntegrator                  (0)
{};

template <class T> NumericalIE<T>::NumericalIE(const NumericalIE& ie)
    : IntegralEngineImplementation<T>(ie)
    , itsMesh         (ie.itsMesh->Clone())
    , itsIntegrator(new MeshIntegrator<T>(itsMesh))
{};

template <class T> NumericalIE<T>::NumericalIE(Mesh* theMesh)
    : IntegralEngineImplementation<T>()
    , itsMesh               (theMesh    )
    , itsIntegrator(new MeshIntegrator<T>(itsMesh))
{};

template <class T> NumericalIE<T>::~NumericalIE()
{
    delete itsIntegrator; //Integrator does not own the mesh.
    delete itsMesh;
}

template <class T> void NumericalIE<T>::Insert(const TBasisSet<T>* theSet)
{
    IntegralEngineImplementation<T>::Insert(theSet);
    itsNormalizations = GetNumericalNormalization(*itsBasisSet);
    assert(!isnan(itsNormalizations));
    CheckInitialized();
}

template <class T> const typename NumericalIE<T>::RVec NumericalIE<T>::GetNumericalNormalization(const TBasisSet<T>& bs) const
{
    CheckInitialized();
    RVec ret=itsIntegrator->Normalize(bs);
    return ret;
}

template <class T> void NumericalIE<T>::CheckInitialized() const
{
    CheckBasisSet();
}

//-----------------------------------------------------------------
//
//  Streamable Object stuff
//
template <class T> std::ostream& NumericalIE<T>::Write(std::ostream& os) const
{
    IntegralEngineImplementation<T>::Write(os);
    return os << itsMesh;
}

template <class T> std::istream& NumericalIE<T>::Read (std::istream& is)
{
    IntegralEngineImplementation<T>::Read(is);
    delete itsMesh;
    itsMesh=Mesh::Factory(is);
    is >> itsMesh;
    delete itsIntegrator;
    itsIntegrator=new MeshIntegrator<T>(itsMesh);
    return is;
}

template <class T> IntegralEngine<T>* NumericalIE<T>::Clone() const
{
    return new NumericalIE(*this);
}

//--------------------------------------------------------------------------------------------
//
//  Overlap type integrals.
//
template <class T> typename NumericalIE<T>::SMat NumericalIE<T>::MakeOverlap() const
{
    CheckInitialized();
    std::cout << "Doing numerical overlap ...";
    std::cout.flush();
    SMat ret=itsIntegrator->Overlap(*itsBasisSet);
    Normalize(ret);
    std::cout << " done" << std::endl;
    return ret;
}

template <class T> typename NumericalIE<T>::Mat NumericalIE<T>::MakeOverlap(const TBasisSet<T>& theOtherBasisSet) const
{
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


template <class T> typename NumericalIE<T>::Vec NumericalIE<T>::MakeOverlap(const ScalarFunction<double>& f) const
{
    CheckInitialized();
    Vec ret=itsIntegrator->Overlap(f,*itsBasisSet);
    assert(!isnan(ret));
    Normalize(ret);
    assert(!isnan(ret));
    return ret;
}


template <class T> void NumericalIE<T>::MakeOverlap3C(MList& mlist, const TBasisSet<T>& bs) const
{
    for (TBasisSetBrowser<T> b(bs); b; b++) mlist.Add(MakeOverlap(*b));
}



//--------------------------------------------------------------------------------------------
//
//  Repulsion type integrals are not implementated, they are 6 dimensional!
//
template <class T> typename NumericalIE<T>::SMat NumericalIE<T>::MakeRepulsion() const
{
    std::cerr << "Error: NumericalIE<T>::MakeRepulsion Do not do repulsion integrals numerically" << std::endl;
    assert(false);
    CheckInitialized();
    SMat ret=itsIntegrator->Repulsion(*itsBasisSet);
    assert(!isnan(ret));
    Normalize(ret);
    assert(!isnan(ret));
    return ret;
}

template <class T> typename NumericalIE<T>::Mat NumericalIE<T>::MakeRepulsion(const TBasisSet<T>& theOtherBasisSet) const
{
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

template <class T> typename NumericalIE<T>::Vec NumericalIE<T>::MakeRepulsion(const ScalarFunction<double>& f) const
{
    CheckInitialized();
    Vec ret=itsIntegrator->Repulsion(f,*itsBasisSet);
    Normalize(ret);
    return ret;
}

template <class T> void NumericalIE<T>::MakeRepulsion3C(MList& ret,const TBasisSet<T>& bs) const
{
    std::cerr << "Error: NumericalIE<T>::MakeRepulsion Do not do repulsion integrals numerically" << std::endl;
    assert(false);
    for (TBasisSetBrowser<T> b(bs); b; b++) ret.Add(MakeRepulsion(*b));
}

template <class T> void NumericalIE<T>::MakeRepulsion4C(ERIList& eris,ERIList&,const BasisGroup*) const
{
    std::cerr << "NumericalIE<T>::MakeRepulsion4C 4 center numerical integrals are not supported" << std::endl;
    exit(-1);
}


//--------------------------------------------------------------------------------
//
//  Special integrals.
//
template <class T> typename NumericalIE<T>::RVec NumericalIE<T>::MakeNormalization() const
{
    CheckInitialized();
    RVec ret=itsNormalizations;
    typename RVec::iterator i(ret.begin());
    BasisSetBrowser           b(*itsBasisSet);
    for(; i!=ret.end()&&b; i++,b++) *i *= b->GetNormalization();
    return ret;
}

template <class T> typename NumericalIE<T>::RVec NumericalIE<T>::MakeCharge() const
{
    CheckInitialized();
    RVec ret=itsIntegrator->Integrate(*itsBasisSet);
    ret=DirectMultiply(ret,itsNormalizations);
    return ret;
}

template <class T> typename NumericalIE<T>::SMat NumericalIE<T>::MakeKinetic() const
{
    CheckInitialized();
    SMat ret=itsIntegrator->Grad(*itsBasisSet);

    Normalize(ret);
    ret  *=  0.5;
    return ret;
}

template <class T> typename NumericalIE<T>::SMat NumericalIE<T>::MakeNuclear(const Cluster& theCluster) const
{
    CheckInitialized();
    SMat ret=itsIntegrator->Nuclear(*itsBasisSet);

    Normalize(ret);
    ret  *=  -(double)(theCluster.GetNuclearCharge());
    return ret;
}

//----------------------------------------------------------
//
//  Private stuff.
//
template <class T> typename NumericalIE<T>::SMat NumericalIE<T>::MakeOverlap(const TBasisFunction<T>& bf) const
{
    CheckInitialized();
    SMat ret=itsIntegrator->Overlap3C(*itsBasisSet,bf);
    Normalize(ret);
    return ret;
}

template <class T> typename NumericalIE<T>::SMat NumericalIE<T>::MakeRepulsion(const TBasisFunction<T>& bf) const
{
    CheckInitialized();
    SMat ret = itsIntegrator->Repulsion3C(*itsBasisSet,bf);
    Normalize(ret);
    return ret;
}


template class NumericalIE<double>;
template class NumericalIE<std::complex<double> >;


