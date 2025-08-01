// File: MeshIntegrator.C  mesh Integrator
module;
#include <iostream>
#include <iomanip>
#include <cmath>
#include <complex>
#include <cassert>

module qchem.Mesh.Integrator;

using std::cout;
using std::endl;

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
template <class T> RVec MeshIntegrator<T>::Integrate(const Vf& v) const
{
    size_t n=v.GetVectorSize();
    RVec ret(n);
    Fill(ret,0.0);

    const Mat& sv(v(*itsMesh));
    int iw=1;
    for (auto rw:*itsMesh)
    {
        for (size_t i=1; i<=n; i++)
            ret(i)+=real(sv(i,iw))*w(rw);
        iw++;
    }

    assert(!isnan(ret));
    return ret;
}

template <class T> RVec MeshIntegrator<T>::Normalize(const Vf& v) const
{
    size_t n=v.GetVectorSize();
    RVec ret(n);
    Fill(ret,0.0);

    const Mat& sv(v(*itsMesh));
    int iw=1;
    for (auto rw:*itsMesh)
    {
        for (size_t i=1; i<=n; i++)
            ret(i)+=real(sv(i,iw)*conj(sv(i,iw))*w(rw));
        iw++;          
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
    size_t n=v.GetVectorSize();
    SMatrix<T> ret(n);
    Fill(ret,T(0.0));

    const Mat& sf(v(*itsMesh));

    int iw=1;
    for (auto rw:*itsMesh)
    {
        for (size_t i=1; i<=n; i++)
            for (size_t j=i; j<=n; j++)
                ret(i,j)+=sf(i,iw)*conj(sf(j,iw))*w(rw);
        iw++;
    }

    return ret;
}

template <class T> typename MeshIntegrator<T>::Vec MeshIntegrator<T>::Overlap(const Rf& f,const Vf& v) const
{
    size_t n=v.GetVectorSize();
    Vec ret(n);
    Fill(ret,T(0.0));

    const Mat & sv(v(*itsMesh));
    assert(!isnan(sv));
    const RVec& sf(f(*itsMesh));
    assert(!isnan(sf));
    int iw=1;
    for (auto rw:*itsMesh)
    {
        for (size_t i=1; i<=n; i++)
            ret(i)+=sv(i,iw)*conj(sf(iw))*w(rw);
        
        iw++;
    }
    assert(!isnan(ret));
    return ret;
}

template <class T> typename MeshIntegrator<T>::Mat MeshIntegrator<T>::Overlap(const Vf& f,const Vf& g) const
{
    size_t nf=f.GetVectorSize();
    size_t ng=g.GetVectorSize();
    Mat ret(nf,ng);
    Fill(ret,T(0.0));

    const Mat& sf(f(*itsMesh));
    const Mat& sg(g(*itsMesh));

    int iw=1;
    for (auto rw:*itsMesh)
    {
        for (size_t fi=1; fi<=nf; fi++)
            for (size_t gi=1; gi<=ng; gi++)
                ret(fi,gi)+=conj(sf(fi,iw))*sg(gi,iw)*w(rw);
        iw++;
    }

    assert(!isnan(ret));
    return ret;
}

//-------------------------------------------------------------------------
//
//  3 center Overlap stuff.
//
template <class T> typename MeshIntegrator<T>::SMat MeshIntegrator<T>::Overlap3C(const Vf& f,const Sf& g) const
{
    size_t n=f.GetVectorSize();
    SMatrix<T> ret(n);
    Fill(ret,T(0.0));

    const Mat& sf(f(*itsMesh));
    const Vec& sg(g(*itsMesh));

    for (size_t i=1; i<=n; i++)
        for (size_t j=i; j<=n; j++)
        {
            int wi=1;
            for (auto rw:*itsMesh)
            {
                ret(i,j)+=conj(sf(i,wi)*sf(j,wi))*sg(wi)*w(rw);
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
template <class T> typename MeshIntegrator<T>::SMat MeshIntegrator<T>::Repulsion(const Vf& f) const
{
    size_t n=f.GetVectorSize();
    SMatrix<T> ret(n);
    Fill(ret,T(0.0));

    const Mat& sf(f(*itsMesh));
    int iw=1;
    for (auto rw1=itsMesh->begin();rw1!=itsMesh->end();rw1++,iw++)
    {
        int jw=iw+1;
        for (auto rw2=rw1+1; rw2!=itsMesh->end(); jw++,rw2++)
        {
            if (r(*rw1)==r(*rw2)) continue;
            double oor=w(*rw1)*w(*rw2)/norm(r(*rw1)-r(*rw2));
            assert(!std::isnan(oor));
            for (size_t i=1; i<=n; i++)
                for (size_t j=i; j<=n; j++)
                    ret(i,j)+=(conj(sf(i,iw))*sf(j,jw)+conj(sf(i,jw))*sf(j,iw))*oor;
        }
    }
    assert(!isnan(ret));
    return ret;
}


template <class T> typename MeshIntegrator<T>::Vec MeshIntegrator<T>::Repulsion(const Rf& h, const Vf& f) const
{
    size_t n=f.GetVectorSize();
    Vec ret(n);
    Fill(ret,T(0.0));

    const Mat & sf(f(*itsMesh));
    const RVec& sh(h(*itsMesh));
    
    int iw=1;
    for (auto rw1=itsMesh->begin();rw1!=itsMesh->end();rw1++,iw++)
    {
        int jw=iw+1;
        for (auto rw2=rw1+1; rw2!=itsMesh->end(); jw++,rw2++)
        {
            if (r(*rw1)==r(*rw2)) continue;
            double oor=w(*rw1)*w(*rw2)/norm(r(*rw1)-r(*rw2));
            assert(!std::isnan(oor));
            for (size_t i=1; i<=n; i++)
                ret(i)+=(sf(i,iw)*conj(sh(jw))+sf(i,jw)*conj(sh(iw)))*oor;
        }
    }
    assert(!isnan(ret));
    return ret;
}


template <class T> typename MeshIntegrator<T>::Mat MeshIntegrator<T>::Repulsion(const Vf& f,const Vf& g) const
{
    size_t nf=f.GetVectorSize();
    size_t ng=g.GetVectorSize();
    Mat ret(nf,ng);
    Fill(ret,T(0.0));

    const Mat& sf(f(*itsMesh));
    const Mat& sg(g(*itsMesh));

    int iw=1;
    for (auto rw1=itsMesh->begin();rw1!=itsMesh->end();rw1++,iw++)
    {
        int jw=iw+1;
        for (auto rw2=rw1+1; rw2!=itsMesh->end(); jw++,rw2++)
        {
            if (r(*rw1)==r(*rw2)) continue;
            double oor=w(*rw1)*w(*rw2)/norm(r(*rw1)-r(*rw2));
            assert(!std::isnan(oor));

            for (size_t i=1; i<=nf; i++)
                for (size_t j=1; j<=ng; j++)
                    ret(i,j)+=(conj(sf(i,iw))*sg(j,jw)+conj(sf(i,jw))*sg(j,iw))*oor;
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
    size_t n=f.GetVectorSize();
    SMatrix<T> ret(n,n);
    Fill(ret,T(0.0));

    const Mat& sf(f(*itsMesh));
    const Vec& sh(h(*itsMesh));

    int iw=1;
    for (auto rw1=itsMesh->begin();rw1!=itsMesh->end();rw1++,iw++)
    {
        int jw=iw+1;
        for (auto rw2=rw1+1; rw2!=itsMesh->end(); jw++,rw2++)
        {
            if (r(*rw1)==r(*rw2)) continue;
            double oor=w(*rw1)*w(*rw2)/norm(r(*rw1)-r(*rw2));
            assert(!std::isnan(oor));
            for (size_t i=1; i<=n; i++)
                for (size_t j=i; j<=n; j++)
                    ret(i,j)+=(conj(sf(i,iw)*sf(j,iw))*sh(jw) + conj(sf(i,jw)*sf(j,jw))*sh(iw))*oor;
        }
    }
    return ret;
}


template <class T> typename MeshIntegrator<T>::SMat MeshIntegrator<T>::Inv_r1(const Vf& f) const
{
    size_t n=f.GetVectorSize();
    SMatrix<T> ret(n,n);
    Fill(ret,T(0.0));

    const Mat& sf(f(*itsMesh));

    for (size_t i=1; i<=n; i++)
        for (size_t j=i; j<=n; j++)
        {
            int wi=1;
            for (auto rw:*itsMesh)
            {
                if (norm(r(rw))!=0) ret(i,j)+=conj(sf(i,wi))*sf(j,wi)*w(rw)/norm(r(rw));
                wi++;
            }
        }
    return ret;
}
template <class T> typename MeshIntegrator<T>::SMat MeshIntegrator<T>::Inv_r2(const Vf& f) const
{
    size_t n=f.GetVectorSize();
    SMatrix<T> ret(n,n);
    Fill(ret,T(0.0));

    const Mat& sf(f(*itsMesh));

    for (size_t i=1; i<=n; i++)
        for (size_t j=i; j<=n; j++)
        {
            int wi=1;
            for (auto rw:*itsMesh)
            {
                double mr=norm(r(rw));
                if (mr!=0) ret(i,j)+=conj(sf(i,wi))*sf(j,wi)*w(rw)/(mr*mr);
                wi++;
            }
        }
    return ret;
}

template <class T> typename MeshIntegrator<T>::SMat MeshIntegrator<T>::Grad(const Vf& f) const
{
    size_t n=f.GetVectorSize();
    SMatrix<T> ret(n);
    Fill(ret,T(0.0));

    const Matrix<Vec3>& sf(f.Gradient(*itsMesh));

    for (size_t i=1; i<=n; i++)
        for (size_t j=i; j<=n; j++)
        {
            int wi=1;
            for (auto rw:*itsMesh)
            {
                 ret(i,j)+=conj(sf(i,wi))*sf(j,wi)*w(rw);
                 wi++;
            }
        }
    return ret;
}

//  Here we use grad(f)=g=g^hat*df/dr, so norm(grad(f))=df/dr which is valid for l=0;
template <class T> typename MeshIntegrator<T>::Mat MeshIntegrator<T>::Grada_b(const Vf& a,const Vf& b) const
{
    size_t na=a.GetVectorSize();
    size_t nb=b.GetVectorSize();
    Mat ret(na,nb);
    Fill(ret,T(0.0));

    const Matrix<Vec3>& sa(a.Gradient(*itsMesh));
    const Mat& sb(b(*itsMesh));

    for (size_t i=1; i<=na; i++)
        for (size_t j=1; j<=nb; j++)
        {
            int wi=1;
            for (auto rw:*itsMesh)
            {
                ret(i,j)+=conj(norm(sa(i,wi)))*sb(j,wi)*w(rw);
                wi++;
            }
        }
    return ret;
}

template <class T> typename MeshIntegrator<T>::Mat MeshIntegrator<T>::a_Gradb(const Vf& a,const Vf& b) const
{
    size_t na=a.GetVectorSize();
    size_t nb=b.GetVectorSize();
    Mat ret(na,nb);
    Fill(ret,T(0.0));

    const Mat& sa(a(*itsMesh));
    const Matrix<Vec3>& sb(b.Gradient(*itsMesh));

    for (size_t i=1; i<=na; i++)
        for (size_t j=1; j<=nb; j++)
        {
            int wi=1;
            for (auto rw:*itsMesh)
            {
                ret(i,j)+=conj(sa(i,wi))*norm(sb(j,wi))*w(rw);
                wi++;
            }
        }
    return ret;
}

template class MeshIntegrator<double>;
template class MeshIntegrator<std::complex<double> >;
