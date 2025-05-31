// File: BSpline/IE_DFT.H  Integrals for BSpline DFT basis sets.
#ifndef _BSpline_IE_DFT_H_
#define _BSpline_IE_DFT_H_

#include "IE.H"
#include "Imp/DataBase/DB_DFT.H"

namespace BSpline
{

template <class T, size_t K> class IE_DFT 
: public virtual Primative_Overlap<T,K>
, public virtual Primative_Repulsion<T,K>
, public DB_DFT<T>
{
    typedef Integrals_Base<T> Base;
    typedef typename Base::fbs_t fbs_t;
    typedef typename Base::SMat SMat;
    typedef typename Base::ERI3 ERI3;
    typedef bspline::Spline<T, K> spline_t;
protected:
    IE_DFT(const DB_cache<T>* db) : DB_DFT<T>(db) {};
    
    virtual ERI3 MakeOverlap3C  (const fbs_t& c) const;
    virtual ERI3 MakeRepulsion3C(const fbs_t& c) const;
private:
    typedef typename BSpline::IrrepIEClient<K>::bf_tuple bf_tuple;
    SMat MakeOverlap  (const bf_tuple& c) const; //ab loops
    SMat MakeRepulsion(const bf_tuple& c) const; //ab loops
};

} //namespace


#endif //_BSpline_IE_DFT_H_
