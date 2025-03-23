// File: NumericalIE.C  Here is where all the numerical integrals get calculated.



#include "BasisSetImplementation/NumericalIEImp.H"
#include "Mesh/MeshIntegrator.H"
#include "BasisSet.H"
#include "Mesh/Mesh.H"
#include "oml/vector.h"
#include <cassert>

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

template <class T> const typename NumericalIEImp<T>::RVec NumericalIEImp<T>::GetNumericalNormalization(bs_t& bs) const
{
    return itsIntegrator->Normalize(bs);
}

template <class T> NumericalIE<T>* NumericalIEImp<T>::Clone() const
{
    return new NumericalIEImp(*this);
}

template <class T> typename NumericalIEImp<T>::Vec NumericalIEImp<T>::MakeOverlap(Vf& bs, Rf& f) const
{
     return itsIntegrator->Overlap(f,bs);
}

template <class T> typename NumericalIEImp<T>::Vec NumericalIEImp<T>::MakeRepulsion(Vf& bs, Rf& f) const
{
    return itsIntegrator->Repulsion(f,bs);
}

template <class T> typename NumericalIEImp<T>::RVec NumericalIEImp<T>::MakeNormalization(bs_t& a) const
{
    return itsIntegrator->Normalize(a);;
}


#ifndef UT_COVERAGE_ONLY


template <class T> typename NumericalIEImp<T>::RVec NumericalIEImp<T>::MakeCharge(bs_t& a) const
{
    //No UT coverage.
    CheckInitialized();
    RVec ret=itsIntegrator->Integrate(a);
    ret=DirectMultiply(ret,itsNormalizations);
    return ret;
}

template <class T> typename NumericalIEImp<T>::SMat NumericalIEImp<T>::MakeOverlap(bs_t& a) const
{
    //No UT coverage.
    CheckInitialized();
    std::cout << "Doing numerical overlap ...";
    std::cout.flush();
    SMat ret=itsIntegrator->Overlap(a);
    Normalize(ret);
    std::cout << " done" << std::endl;
    return ret;
}

template <class T> typename NumericalIEImp<T>::Mat NumericalIEImp<T>::MakeOverlap(bs_t& a,bs_t& b) const
{
    //No UT coverage.
    CheckInitialized();
//
//  Can't assume other basis set is numerical, so we must explicitly evaluate
//  the normalization constants over this mesh.
//
    RVec otherNormalizations = GetNumericalNormalization(b);
    Mat ret=itsIntegrator->Overlap(a,b);
    Normalize(itsNormalizations,ret,otherNormalizations);
    return ret;
}


template <class T> typename NumericalIEImp<T>::Mat NumericalIEImp<T>::MakeRepulsion(bs_t& a,bs_t& b    ) const
{
    //No UT coverage.
    std::cerr << "Error: NumericalIE<T>::MakeRepulsion Do not do repulsion integrals numerically" << std::endl;
    assert(false);
    CheckInitialized();
//
//  Can't assume other basis set is numerical, so we must explicitly evaluate
//  the normalization constants over this mesh.
//
    RVec otherNormalizations = GetNumericalNormalization(b);
    Mat ret=itsIntegrator->Repulsion(a,b);
    Normalize(itsNormalizations,ret,otherNormalizations);
    return ret;
}


template <class T> typename NumericalIEImp<T>::SMat NumericalIEImp<T>::MakeRepulsion(bs_t& a ) const
{
    //No UT coverage.
    std::cerr << "Error: NumericalIE<T>::MakeRepulsion Do not do repulsion integrals numerically" << std::endl;
    assert(false);
    CheckInitialized();
    SMat ret=itsIntegrator->Repulsion(a);
    assert(!isnan(ret));
    Normalize(ret);
    assert(!isnan(ret));
    return ret;
}
#endif

#ifdef USE_FOR_DEBUGGING_ANALYTIC

template <class T> void NumericalIEImp<T>::MakeOverlap3C(ERI3& mlist, bs_t& bs) const
{
    // No UT coverage
    for (auto b=bs.beginT(); b!=bs.end(); b++) mlist.Add(MakeOverlap(**b));
}

template <class T> void NumericalIEImp<T>::MakeRepulsion3C(ERI3& ret,bs_t& bs) const
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
template <class T> typename NumericalIEImp<T>::SMat NumericalIEImp<T>::MakeGrad2() const
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


template class NumericalIEImp<double>;
template class NumericalIEImp<std::complex<double> >;


