// File: BSpline/IE.H Common IE code for BSpline basis sets.
#ifndef _BSpline_IE_H_
#define _BSpline_IE_H_

#include "../../BFGrouper.H"
#include "Imp/Integrals/Cache4.H"
#include "Imp/DataBase/HeapDB.H"
#include <bspline/Core.h>


namespace BSpline
{

// K is the spline order
template <class T, size_t K> class Primative_Overlap
{
    typedef bspline::Spline<T, K> spline_t;
public:
    virtual double Overlap(const spline_t& a , const spline_t& b,size_t l_total) const=0;
};
template <class T, size_t K> class Primative_Grad2
{
    typedef bspline::Spline<T, K> spline_t;
public:
    virtual double Grad2(const spline_t& a , const spline_t& b,size_t la, size_t lb) const=0;
};
template <class T, size_t K> class Primative_Nuclear
{
    typedef bspline::Spline<T, K> spline_t;
public:
    virtual double Nuclear(const spline_t& a , const spline_t& b,size_t l_total) const=0;
};
template <class T, size_t K> class Primative_Repulsion
{
    typedef bspline::Spline<T, K> spline_t;
public:
    virtual double Repulsion(const spline_t& a , const spline_t& b,size_t la, size_t lc) const=0;
};
template <class T, size_t K> class Primative_Charge
{
    typedef bspline::Spline<T, K> spline_t;
public:
    virtual double Charge   (const spline_t& a, size_t l) const=0;
};

template <class T, size_t K> class IE_Overlap
: public virtual Primative_Overlap<T,K>
, public DB_Overlap<T>
{
protected:
    typedef bspline::Spline<T, K> spline_t;
    using Primative_Overlap<T,K>::Overlap;
    virtual typename Integrals_Base<T>::SMat MakeOverlap() const;
    IE_Overlap(const DB_cache<T>* db) : DB_Overlap<T>(db) {};
};
template <class T, size_t K> class IE_Grad2
: public virtual Primative_Grad2<T,K>
, public DB_Kinetic<T>
{
protected:
    using Primative_Grad2<T,K>::Grad2;
    virtual typename Integrals_Base<T>::SMat MakeKinetic() const;
    IE_Grad2(const DB_cache<T>* db) : DB_Kinetic<T>(db) {};
};
template <class T, size_t K> class IE_Nuclear
: public virtual Primative_Nuclear<T,K>
, public DB_Nuclear<T>
{
protected:
    using Primative_Nuclear<T,K>::Nuclear;
    virtual typename Integrals_Base<T>::SMat MakeNuclear(const Cluster* cl) const;
    IE_Nuclear(const DB_cache<T>* db) : DB_Nuclear<T>(db) {};
};
template <class T, size_t K> class IE_XGrad2
: public virtual Primative_Grad2<T,K>
, public DB_XGrad2<T>
{
protected:
    using Primative_Grad2<T,K>::Grad2;
    virtual typename Integrals_Base<T>::Mat MakeKinetic(const Orbital_RKBS_IBS<T>* rkbs) const;
    IE_XGrad2(const DB_cache<T>* db) : DB_XGrad2<T>(db) {};
};

} //namespace
#endif // _BSpline_IE_H_
