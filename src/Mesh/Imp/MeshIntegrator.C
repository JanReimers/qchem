// File: MeshIntegrator.C  mesh Integrator
module;
#include <iostream>
#include <iomanip>
#include <cmath>
#include <complex>
#include <cassert>
#include "blaze/Math.h"
module qchem.Mesh.Integrator;
import qchem.Conversions;
import qchem.Vector3D;
using std::cout;
using std::endl;

namespace std
{
    double conj(double a) {return a;}
    double real(double a) {return a;}
}
//-------------------------------------------------------------------------
//
//  Construction zone.
//
template <class T> MeshIntegrator<T>::MeshIntegrator(const Mesh* m)
    : itsMesh(m)
{
    assert(itsMesh);
};

//-------------------------------------------------------------------------
//
//  Straight integration.
//
template <class T> rvec_t MeshIntegrator<T>::Integrate(const Vf& v) const
{
    size_t n=v.GetVectorSize();
    rvec_t ret(n,0.0);

    const mat_t<T>& sv(v(*itsMesh));
    int iw=0;
    for (auto rw:*itsMesh)
    {
        for (size_t i=0; i<n; i++)
            ret[i]+=std::real(sv(i,iw))*w(rw);
        iw++;
    }

    assert(!isnan(ret));
    return ret;
}

template <class T> rvec_t MeshIntegrator<T>::Normalize(const Vf& v) const
{
    size_t n=v.GetVectorSize();
    rvec_t ret(n,0.0);

    const mat_t<T>& sv(v(*itsMesh));
    int iw=0;
    for (auto rw:*itsMesh)
    {
        for (size_t i=0; i<n; i++)
            ret[i]+=std::real(sv(i,iw)*std::conj(sv(i,iw))*w(rw));
        iw++;          
    }

    ret=1.0/sqrt(ret);
    return ret;
}

//-------------------------------------------------------------------------
//
//  2 center Overlap stuff.
//
template <class T> smat_t<T> MeshIntegrator<T>::Overlap(const Vf& v) const
{
    size_t n=v.GetVectorSize();
    smat_t<T> ret=zero<T>(n);

    const mat_t<T>& sf(v(*itsMesh));
    int iw=0;
    for (auto rw:*itsMesh)
    {
        for (size_t i=0; i<n; i++)
            for (size_t j=i; j<n; j++)
                ret(i,j)+=sf(i,iw)*std::conj(sf(j,iw))*w(rw);
        iw++;
    }

    return ret;
}

template <class T> vec_t<T> MeshIntegrator<T>::Overlap(const Rf& f,const Vf& v) const
{
    size_t n=v.GetVectorSize();
    vec_t<T> ret(n,T(0.0));

    const mat_t<T> & sv(v(*itsMesh));
    assert(!isnan(sv));
    const rvec_t& sf(f(*itsMesh));
    assert(!isnan(sf));
    int iw=0;
    for (auto rw:*itsMesh)
    {
        for (size_t i=0; i<n; i++)
            ret[i]+=sv(i,iw)*std::conj(sf[iw])*w(rw);
        
        iw++;
    }
    assert(!isnan(ret));
    return ret;
}

template <class T> mat_t<T> MeshIntegrator<T>::Overlap(const Vf& f,const Vf& g) const
{
    size_t nf=f.GetVectorSize();
    size_t ng=g.GetVectorSize();
    mat_t<T> ret(nf,ng,T(0.0));

    const mat_t<T>& sf(f(*itsMesh));
    const mat_t<T>& sg(g(*itsMesh));

    int iw=0;
    for (auto rw:*itsMesh)
    {
        for (size_t fi=0; fi<nf; fi++)
            for (size_t gi=0; gi<ng; gi++)
                ret(fi,gi)+=std::conj(sf(fi,iw))*sg(gi,iw)*w(rw);
        iw++;
    }

    assert(!isnan(ret));
    return ret;
}

//-------------------------------------------------------------------------
//
//  3 center Overlap stuff.
//
template <class T> smat_t<T> MeshIntegrator<T>::Overlap3C(const Vf& f,const Sf& g) const
{
    size_t n=f.GetVectorSize();
    smat_t<T> ret=zero<T>(n);

    const mat_t<T>& sf(f(*itsMesh));
    const vec_t<T>& sg(g(*itsMesh));

    for (size_t i=0; i<n; i++)
        for (size_t j=0; j<n; j++)
        {
            int wi=0;
            for (auto rw:*itsMesh)
            {
                ret(i,j)+=std::conj(sf(i,wi)*sf(j,wi))*sg[wi]*w(rw);
                wi++;
            }
        }

    assert(!isnan(ret));
    return ret;
}

//---------------------------------------------------------------------------------
//
//  2 center Repulsion stuff.  These are usually not very accurate.
//
template <class T> smat_t<T> MeshIntegrator<T>::Repulsion(const Vf& f) const
{
    size_t n=f.GetVectorSize();
    smat_t<T> ret=zero<T>(n);

    const mat_t<T>& sf(f(*itsMesh));
    int iw=0;
    for (auto rw1=itsMesh->begin();rw1!=itsMesh->end();rw1++,iw++)
    {
        int jw=iw+1;
        for (auto rw2=rw1+1; rw2!=itsMesh->end(); jw++,rw2++)
        {
            if (r(*rw1)==r(*rw2)) continue;
            double oor=w(*rw1)*w(*rw2)/norm(r(*rw1)-r(*rw2));
            assert(!std::isnan(oor));
            for (size_t i=0; i<n; i++)
                for (size_t j=i; j<n; j++)
                    ret(i,j)+=(std::conj(sf(i,iw))*sf(j,jw)+std::conj(sf(i,jw))*sf(j,iw))*oor;
        }
    }
    assert(!isnan(ret));
    return ret;
}


template <class T> vec_t<T> MeshIntegrator<T>::Repulsion(const Rf& h, const Vf& f) const
{
    size_t n=f.GetVectorSize();
    vec_t<T> ret(n,T(0.0));

    const mat_t<T> & sf(f(*itsMesh));
    const rvec_t& sh(h(*itsMesh));
    
    int iw=0;
    for (auto rw1=itsMesh->begin();rw1!=itsMesh->end();rw1++,iw++)
    {
        int jw=iw+1;
        for (auto rw2=rw1+1; rw2!=itsMesh->end(); jw++,rw2++)
        {
            if (r(*rw1)==r(*rw2)) continue;
            double oor=w(*rw1)*w(*rw2)/norm(r(*rw1)-r(*rw2));
            assert(!std::isnan(oor));
            for (size_t i=0; i<n; i++)
                ret[i]+=(sf(i,iw)*std::conj(sh[jw])+sf(i,jw)*std::conj(sh[iw]))*oor;
        }
    }
    assert(!isnan(ret));
    return ret;
}


template <class T> mat_t<T> MeshIntegrator<T>::Repulsion(const Vf& f,const Vf& g) const
{
    size_t nf=f.GetVectorSize();
    size_t ng=g.GetVectorSize();
    mat_t<T> ret(nf,ng,T(0.0));

    const mat_t<T>& sf(f(*itsMesh));
    const mat_t<T>& sg(g(*itsMesh));

    int iw=0;
    for (auto rw1=itsMesh->begin();rw1!=itsMesh->end();rw1++,iw++)
    {
        int jw=iw+1;
        for (auto rw2=rw1+1; rw2!=itsMesh->end(); jw++,rw2++)
        {
            if (r(*rw1)==r(*rw2)) continue;
            double oor=w(*rw1)*w(*rw2)/norm(r(*rw1)-r(*rw2));
            assert(!std::isnan(oor));

            for (size_t i=0; i<nf; i++)
                for (size_t j=0; j<ng; j++)
                    ret(i,j)+=(std::conj(sf(i,iw))*sg(j,jw)+std::conj(sf(i,jw))*sg(j,iw))*oor;
        }
    }
    assert(!isnan(ret));
    return ret;
}


//---------------------------------------------------------------------------------
//
//  3 center Repulsion stuff.  These are usually not very accurate.
//
template <class T> smat_t<T> MeshIntegrator<T>::Repulsion3C(const Vf& f, const Sf& h) const
{
    size_t n=f.GetVectorSize();
    smat_t<T> ret=zero<T>(n);

    const mat_t<T>& sf(f(*itsMesh));
    const vec_t<T>& sh(h(*itsMesh));

    int iw=0;
    for (auto rw1=itsMesh->begin();rw1!=itsMesh->end();rw1++,iw++)
    {
        int jw=iw+1;
        for (auto rw2=rw1+1; rw2!=itsMesh->end(); jw++,rw2++)
        {
            if (r(*rw1)==r(*rw2)) continue;
            double oor=w(*rw1)*w(*rw2)/norm(r(*rw1)-r(*rw2));
            assert(!std::isnan(oor));
            for (size_t i=0; i<n; i++)
                for (size_t j=0; j<n; j++)
                    ret(i,j)+=(std::conj(sf(i,iw)*sf(j,iw))*sh[jw] + std::conj(sf(i,jw)*sf(j,jw))*sh[iw])*oor;
        }
    }
    return ret;
}


template <class T> smat_t<T> MeshIntegrator<T>::Inv_r1(const Vf& f) const
{
    size_t n=f.GetVectorSize();
    smat_t<T> ret=zero<T>(n);

    const mat_t<T>& sf(f(*itsMesh));

    for (size_t i=0; i<n; i++)
        for (size_t j=i; j<n; j++)
        {
            int wi=0;
            for (auto rw:*itsMesh)
            {
                if (norm(r(rw))!=0) ret(i,j)+=std::conj(sf(i,wi))*sf(j,wi)*w(rw)/norm(r(rw));
                wi++;
            }
        }
    return ret;
}
template <class T> smat_t<T> MeshIntegrator<T>::Inv_r2(const Vf& f) const
{
    size_t n=f.GetVectorSize();
    smat_t<T> ret=zero<T>(n);

    const mat_t<T>& sf(f(*itsMesh));

    for (size_t i=0; i<n; i++)
        for (size_t j=i; j<n; j++)
        {
            int wi=0;
            for (auto rw:*itsMesh)
            {
                double mr=norm(r(rw));
                if (mr!=0) ret(i,j)+=std::conj(sf(i,wi))*sf(j,wi)*w(rw)/(mr*mr);
                wi++;
            }
        }
    return ret;
}

template <class T> smat_t<T> MeshIntegrator<T>::Grad2(const Vf& f) const
{
    size_t n=f.GetVectorSize();
    smat_t<T> ret=zero<T>(n);

    const mat_t<Vec3>& sf(f.Gradient(*itsMesh));

    for (size_t i=0; i<n; i++)
        for (size_t j=i; j<n; j++)
        {
            int wi=0;
            for (auto rw:*itsMesh)
            {
                 ret(i,j)+=conj(sf(i,wi))*sf(j,wi)*w(rw);
                 wi++;
            }
        }
    return ret;
}

//  Here we use grad(f)=g=g^hat*df/dr, so norm(grad(f))=df/dr which is valid for l=0;
template <class T> mat_t<T> MeshIntegrator<T>::Grada_b(const Vf& a,const Vf& b) const
{
    size_t na=a.GetVectorSize();
    size_t nb=b.GetVectorSize();
    mat_t<T> ret(na,nb,T(0.0));
 
    const mat_t<Vec3>& sa(a.Gradient(*itsMesh));
    const mat_t<T>& sb(b(*itsMesh));

    for (size_t i=0; i<na; i++)
        for (size_t j=0; j<nb; j++)
        {
            int wi=0;
            for (auto rw:*itsMesh)
            {
                ret(i,j)+=std::conj(norm(sa(i,wi)))*sb(j,wi)*w(rw);
                wi++;
            }
        }
    return ret;
}

template <class T> mat_t<T> MeshIntegrator<T>::a_Gradb(const Vf& a,const Vf& b) const
{
    size_t na=a.GetVectorSize();
    size_t nb=b.GetVectorSize();
    mat_t<T> ret(na,nb,T(0.0));
   
    const mat_t<T>& sa(a(*itsMesh));
    const mat_t<Vec3>& sb(b.Gradient(*itsMesh));

    for (size_t i=0; i<na; i++)
        for (size_t j=0; j<nb; j++)
        {
            int wi=0;
            for (auto rw:*itsMesh)
            {
                ret(i,j)+=std::conj(sa(i,wi))*norm(sb(j,wi))*w(rw);
                wi++;
            }
        }
    return ret;
}

template class MeshIntegrator<double>;
template class MeshIntegrator<std::complex<double> >;
