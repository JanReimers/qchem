// File: BSpline/IE_DHF.H  Integrals for BSpline DHF basis sets.
#ifndef _BSpline_IE_DHF_H_
#define _BSpline_IE_DHF_H_

#include "IE.H"
#include "DataBase/DB_DHF.H"

namespace BSpline
{

template <class T, size_t K> class IE_RKBL 
    : public IE_Overlap<T,K>
    , public IE_XGrad2 <T,K>
    , public IE_Inv_r1<T,K>
{
protected:
    IE_RKBL(const DB_cache<T>* db) : IE_Overlap<T,K>(db),IE_XGrad2<T,K>(db),IE_Inv_r1<T,K>(db) {};

};
template <class T, size_t K> class IE_RKBS 
: public IE_Kinetic  <T,K>
, public IE_Inv_r1<T,K>
{
protected:
    IE_RKBS(const DB_cache<T>* db) : IE_Kinetic<T,K>(db), IE_Inv_r1<T,K>(db) {};
};

} //namespace

#endif //_BSpline_IE_DHF_H_
