// File: BSpline/IEC.H Common IE client code for all atom BSpline basis sets and IEs.
#ifndef _BSplineIEClient_H_
#define _BSplineIEClient_H_

#include <IEClient.H> // <> prevents file search on this folder.
#include <bspline/Core.h>
#include "oml/vector.h"

namespace BSpline
{
template <size_t K> struct IrrepIEClient : public virtual ::IrrepIEClient
{
    typedef bspline::Spline<double, K> spline_t;
    std::vector<spline_t> splines;
    Vector<double> ns; //Normalization constants
    std::vector<size_t> sp_indices; //Unique spline index

    double rmin,rmax; //This might be needed for creating fit basis sets.
    size_t l;
    int m;
    
    IrrepIEClient(size_t Ngrid, double rmin, double rmax, size_t l, int m=0);
    virtual ~IrrepIEClient();
    
    virtual size_t size() const {return splines.size();}
    auto indices() const {return ns.indices();}
    auto indices(size_t i) const {return ns.indices(i);}
    const spline_t& operator()(index_t i) const {assert(i>0);return splines[i-1];}

    typedef std::tuple<int,int,const spline_t*,double> bf_tuple;
    bf_tuple tuple(index_t i) const {return std::make_tuple(l,m,&splines[i-1],ns(i));}
private:
    IrrepIEClient(const IrrepIEClient<K>&) =delete;
    IrrepIEClient<K>& operator=(const IrrepIEClient<K>&) =delete;
    std::vector<double> MakeLogKnots(size_t NGrid, double rmin, double rmax);
};
} //namespace

#endif // _BSplineIEClient_H_
