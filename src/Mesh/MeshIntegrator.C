// File: MeshIntegrator.C  mesh Integrator


#include "Mesh/MeshIntegrator.H"
#include "Mesh/Mesh.H"
#include "Mesh/MeshImplementation.H"
#include "Mesh/MeshBrowser.H"
#include "Functions/ScalarFunction.H"
#include "Functions/VectorFunction.H"
#include "oml/matrix.h"
#include "oml/smatrix.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <cassert>

//-------------------------------------------------------------------------
//
//  Construction zone.
//
template <class T> MeshIntegrator<T>::MeshIntegrator(const Mesh* theMesh)
    : itsMesh(dynamic_cast<const MeshImplementation*>(theMesh))
{
    assert(itsMesh);
};

//-------------------------------------------------------------------------
//
//  Straight integration.
//
template <class T> typename MeshIntegrator<T>::RVec MeshIntegrator<T>::Integrate(const Vf& v) const
{
    index_t n=v.GetVectorSize();
    RVec ret(n);
    Fill(ret,0.0);

    typename RVec::Subscriptor      sr(ret);
    const Mat& sv(v(*itsMesh));
    MeshBrowser m(*itsMesh);
    for (index_t iw=1; m; iw++,m++)
        for (index_t i=1; i<=n; i++)
            sr(i)+=sv(i,iw)*m.W();

    assert(!isnan(ret));
    return ret;
}

template <class T> typename MeshIntegrator<T>::RVec MeshIntegrator<T>::Normalize(const Vf& v) const
{
    index_t n=v.GetVectorSize();
    RVec ret(n);
    Fill(ret,0.0);

//    v(*itsMesh);
    typename RVec::Subscriptor      sr(ret);
    const Mat& sv(v(*itsMesh));
    MeshBrowser m(*itsMesh);
    for (index_t iw=1; m; iw++,m++)
        for (index_t i=1; i<=n; i++)
        {
//            std::cout << sv(i,iw) << " " << m.W() << std::endl;
            sr(i)+=sv(i,iw)*sv(i,iw)*m.W();

        }

    ret=1.0/sqrt(ret);
    return ret;
}

//-------------------------------------------------------------------------
//
//  2 center Overlap stuff.
//
template <class T> typename MeshIntegrator<T>::SMat MeshIntegrator<T>::Overlap(const Vf& v) const
{
    index_t n=v.GetVectorSize();
    SMat ret(n,n);
    Fill(ret,0.0);

    typename SMat::Subscriptor      r (ret);
    const Mat& sf(v(*itsMesh));

    MeshBrowser m(*itsMesh);
    for (index_t wi=1; m; wi++,m++)
        for (index_t i=1; i<=n; i++)
            for (index_t j=i; j<=n; j++)
                r(i,j)+=sf(i,wi)*sf(j,wi)*m.W();

    return ret;
}

template <class T> typename MeshIntegrator<T>::Vec MeshIntegrator<T>::Overlap(const Rf& f,const Vf& v) const
{
    index_t n=v.GetVectorSize();
    Vec ret(n);
    Fill(ret,0.0);

    typename Vec ::Subscriptor      sr(ret);
    const Mat & sv(v(*itsMesh));
    assert(!isnan(sv));
    const RVec& sf(f(*itsMesh));
    assert(!isnan(sf));
    MeshBrowser m(*itsMesh);
    for (index_t iw=1; m; iw++,m++)
        for (index_t i=1; i<=n; i++)
        {
            sr(i)+=sv(i,iw)*sf(iw)*m.W();
        }
    assert(!isnan(ret));
    return ret;

//  return v(*itsMesh) * conj(DirectMultiply(f(*itsMesh),itsMesh->itsWeights));
}

template <class T> typename MeshIntegrator<T>::Mat MeshIntegrator<T>::Overlap(const Vf& f,const Vf& g) const
{
    index_t nf=f.GetVectorSize();
    index_t ng=g.GetVectorSize();
    Mat ret(nf,ng);
    Fill(ret,0.0);

    typename Mat::Subscriptor      r (ret);
    const Mat& sf(f(*itsMesh));
    const Mat& sg(g(*itsMesh));

    MeshBrowser m(*itsMesh);
    for (index_t wi=1; m; wi++,m++)
        for (index_t fi=1; fi<=nf; fi++)
            for (index_t gi=1; gi<=ng; gi++)
                r(fi,gi)+=sf(fi,wi)*sg(gi,wi)*m.W();

    assert(!isnan(ret));
    return ret;
}

//-------------------------------------------------------------------------
//
//  3 center Overlap stuff.
//
template <class T> typename MeshIntegrator<T>::SMat MeshIntegrator<T>::Overlap3C(const Vf& f,const Sf& g) const
{
    index_t n=f.GetVectorSize();
    SMat ret(n,n);
    Fill(ret,0.0);

    typename SMat::Subscriptor      r (ret);
    const Mat& sf(f(*itsMesh));
    const Vec& sg(g(*itsMesh));

    for (index_t i=1; i<=n; i++)
        for (index_t j=i; j<=n; j++)
        {
            MeshBrowser m(*itsMesh);
            for (index_t wi=1; m; wi++,m++) r(i,j)+=sf(i,wi)*sf(j,wi)*sg(wi)*m.W();
        }

    assert(!isnan(ret));
    return ret;
}

//---------------------------------------------------------------------------------
//
//  2 center Repulsion stuff.  These are usually not very accurate.
//
template <class T> typename MeshIntegrator<T>::SMat MeshIntegrator<T>::Repulsion(const Vf& f) const
{
    index_t n=f.GetVectorSize();
    SMat ret(n,n);
    Fill(ret,0.0);

    typename SMat::Subscriptor          r (ret);
    const Mat& sf(f(*itsMesh));

    MeshBrowser m1(*itsMesh);
//    double oor1=m1.W()*m1.W()/!m1.R();
    for (index_t wi=1; m1; wi++,m1++)
    {
//        for (index_t i=1; i<=n; i++)
//            for (index_t j=i; j<=n; j++)
//                r(i,j)+=sf(i,wi)*sf(j,wi)*oor1;

        MeshBrowser m2(m1);
        m2++;
        double oor=m1.W()*m2.W()/!(m1.R()-m2.R());
        assert(!std::isnan(oor));
        for (index_t wj=wi+1; m2; wj++,m2++)
        {
            for (index_t i=1; i<=n; i++)
                for (index_t j=i; j<=n; j++)
                    r(i,j)+=(sf(i,wi)*sf(j,wj)+sf(i,wj)*sf(j,wi))*oor;
        }
    }
    assert(!isnan(ret));
    return ret;
}


template <class T> typename MeshIntegrator<T>::Vec MeshIntegrator<T>::Repulsion(const Rf& h, const Vf& f) const
{
    index_t n=f.GetVectorSize();
    Vec ret(n);
    Fill(ret,0.0);

    typename Vec ::Subscriptor      r (ret);
    const Mat & sf(f(*itsMesh));
    const RVec& sh(h(*itsMesh));

    MeshBrowser m1(*itsMesh);
    for (index_t wi=1; m1; wi++,m1++)
    {
        double oor1=m1.W()*m1.W()/!m1.R();
        for (index_t i=1; i<=n; i++)
            r(i)+=sf(i,wi)*sh(wi)*oor1;

        MeshBrowser m2(m1);
        m2++;
        for (index_t wj=wi+1; m2; wj++,m2++)
        {
            double oor=m1.W()*m2.W()/!(m1.R()-m2.R());
            for (index_t i=1; i<=n; i++)
                r(i)+=(sf(i,wi)*sh(wj) + sf(i,wj)*sh(wi))*oor;
        }
    }
    assert(!isnan(ret));
    return ret;
}

template <class T> typename MeshIntegrator<T>::Mat MeshIntegrator<T>::Repulsion(const Vf& f,const Vf& g) const
{
    index_t nf=f.GetVectorSize();
    index_t ng=g.GetVectorSize();
    Mat ret(nf,ng);
    Fill(ret,0.0);

    typename Mat::Subscriptor      r (ret);
    const Mat& sf(f(*itsMesh));
    const Mat& sg(g(*itsMesh));

    MeshBrowser m1(*itsMesh);
    for (index_t wi=1; m1; wi++,m1++)
    {
        double oor1=m1.W()*m1.W()/!m1.R();
        for (index_t i=1; i<=nf; i++)
            for (index_t j=1; j<=ng; j++)
                r(i,j)+=sf(i,wi)*sg(j,wi)*oor1;

        MeshBrowser m2(m1);
        m2++;
        for (index_t wj=wi+1; m2; wj++,m2++)
        {
            double oor=m1.W()*m2.W()/!(m1.R()-m2.R());
            for (index_t i=1; i<=nf; i++)
                for (index_t j=1; j<=ng; j++)
                    r(i,j)+=(sf(i,wi)*sg(j,wj)+sf(i,wj)*sg(j,wi))*oor;
        }
    }
    assert(!isnan(ret));
    return ret;
}


//---------------------------------------------------------------------------------
//
//  3 center Repulsion stuff.  These are usually not very accurate.
//
template <class T> typename MeshIntegrator<T>::SMat MeshIntegrator<T>::Repulsion3C(const Vf& f, const Sf& h) const
{
    index_t n=f.GetVectorSize();
    SMat ret(n,n);
    Fill(ret,0.0);

    typename SMat::Subscriptor      r (ret);
    const Mat& sf(f(*itsMesh));
    const Vec& sh(h(*itsMesh));

    MeshBrowser m1(*itsMesh);
    for (index_t wi=1; m1; wi++,m1++)
    {
        double oor1=m1.W()*m1.W()/!m1.R();
        for (index_t i=1; i<=n; i++)
            for (index_t j=i; j<=n; j++)
                r(i,j)+=sf(i,wi)*sf(j,wi)*sh(wi)*oor1;

        MeshBrowser m2(m1);
        m2++;
        for (index_t wj=wi+1; m2; wj++,m2++)
        {
            double oor=m1.W()*m2.W()/!(m1.R()-m2.R());
            for (index_t i=1; i<=n; i++)
                for (index_t j=i; j<=n; j++)
                    r(i,j)+=(sf(i,wi)*sf(j,wi)*sh(wj) + sf(i,wj)*sf(j,wj)*sh(wi))*oor;
        }
    }
    return ret;
}


template <class T> typename MeshIntegrator<T>::SMat MeshIntegrator<T>::Nuclear(const Vf& f) const
{
    index_t n=f.GetVectorSize();
    SMat ret(n,n);
    Fill(ret,0.0);

    typename SMat::Subscriptor   r (ret);
    const Mat& sf(f(*itsMesh));

    for (index_t i=1; i<=n; i++)
        for (index_t j=i; j<=n; j++)
        {
            MeshBrowser m(*itsMesh);
            for (index_t wi=1; m; wi++,m++) if (!m.R()!=0) r(i,j)+=sf(i,wi)*sf(j,wi)*m.W()/!m.R();
        }
    return ret;
}

template <class T> typename MeshIntegrator<T>::SMat MeshIntegrator<T>::Grad(const Vf& f) const
{
    index_t n=f.GetVectorSize();
    SMat ret(n,n);
    Fill(ret,0.0);

    typename SMat::Subscriptor      r (ret);
    const Matrix<Vec3>& sf(f.Gradient(*itsMesh));

    for (index_t i=1; i<=n; i++)
        for (index_t j=i; j<=n; j++)
        {
            MeshBrowser m(*itsMesh);
            for (index_t wi=1; m; wi++,m++) r(i,j)+=sf(i,wi)*sf(j,wi)*m.W();
        }
    return ret;
}

template class MeshIntegrator<double>;
