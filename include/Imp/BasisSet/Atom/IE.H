// File: AtomIE.H Common IE code for all atom basis sets.
#ifndef _AtomIE_H_
#define _AtomIE_H_

#include "BFGrouper.H"
#include "Imp/Integrals/Cache4.H"
#include "Imp/DataBase/HeapDB.H"


template <class T> class Primative_Overlap
{
public:
    virtual double Overlap(double ea ,double eb,size_t l_total) const=0;
};
template <class T> class Primative_Grad2
{
public:
    virtual double Grad2(double ea ,double eb,size_t la, size_t lb) const=0;
};
template <class T> class Primative_Inv_r1
{
public:
    virtual double Inv_r1(double ea ,double eb,size_t l_total) const=0;
};
template <class T> class Primative_Inv_r2
{
public:
    virtual double Inv_r2(double ea ,double eb,size_t l_total) const=0;
};
template <class T> class Primative_Repulsion
{
public:
    virtual double Repulsion(double ea ,double ec,size_t la, size_t lc) const=0;
};
template <class T> class Primative_Charge
{
public:
    virtual double Charge   (double ea, size_t l) const=0;
};

template <class T> class AtomIE_Overlap
: public virtual Primative_Overlap<T>
, public DB_Overlap<T>
{
protected:
    using Primative_Overlap<T>::Overlap;
    virtual typename Integrals_Base<T>::SMat MakeOverlap() const;
    AtomIE_Overlap(const DB_cache<T>* db) : DB_Overlap<T>(db) {};
};
template <class T> class AtomIE_Kinetic
: public virtual Primative_Grad2<T>
, public virtual Primative_Inv_r2<T> //for centrifugal term.
, public DB_Kinetic<T>
{
protected:
    using Primative_Grad2<T>::Grad2;
    virtual typename Integrals_Base<T>::SMat MakeKinetic() const;
    AtomIE_Kinetic(const DB_cache<T>* db) : DB_Kinetic<T>(db) {};
};
template <class T> class AtomIE_Nuclear
: public virtual Primative_Inv_r1<T>
, public DB_Nuclear<T>
{
protected:
    using Primative_Inv_r1<T>::Inv_r1;
    virtual typename Integrals_Base<T>::SMat MakeNuclear(const Cluster* cl) const;
    AtomIE_Nuclear(const DB_cache<T>* db) : DB_Nuclear<T>(db) {};
};
template <class T> class AtomIE_XGrad2
: public virtual Primative_Grad2<T>
, public DB_XGrad2<T>
{
protected:
    using Primative_Grad2<T>::Grad2;
    virtual typename Integrals_Base<T>::Mat MakeKinetic(const Orbital_RKBS_IBS<T>* rkbs) const;
    AtomIE_XGrad2(const DB_cache<T>* db) : DB_XGrad2<T>(db) {};
};

#endif // _AtomIE_H_
