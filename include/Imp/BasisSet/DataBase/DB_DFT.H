// File: DB_DFT.H Cache DFT 3 center integrals.
#ifndef _DB_DFT_H_
#define _DB_DFT_H_

#include "HeapDB.H"
#include <DFT_IE.H>

template <class T> class DB_DFT 
    : virtual public Integrals_DFT<T>
    , public DB_Common<T>
{
    typedef typename Integrals_Base<T>::fbs_t fbs_t;
    typedef typename Integrals_Base<T>::ERI3 ERI3;
protected:
    DB_DFT(const DB_cache<T>* db) : DB_Common<T>(db) {}
    virtual const ERI3& Overlap3C  (const fbs_t& c) const; //<ab|c>
    virtual const ERI3& Repulsion3C(const fbs_t& c) const; //<a(1)b(1)|1/r12|c(2)>
    virtual ERI3 MakeOverlap3C  (const fbs_t& c) const=0;
    virtual ERI3 MakeRepulsion3C(const fbs_t& c) const=0;
};

#endif //_DB_DFT_H_