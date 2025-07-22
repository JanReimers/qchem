// File: BSpline/IEC.C Common IE client code for all atom BSpline basis sets and IEs.
module;
#include <bspline/Core.h>
#include "radial/BSpline/GLQuadrature.H"
export module qchem.Basisset.Atom.radial.BSpline.IEC;
import qchem.BasisSet.Atom.IEClient;
import oml;

export namespace BSpline
{
template <size_t K> struct IrrepIEClient 
    : public virtual ::IrrepIEClient
    , public AtomIrrepIEClient
{
    typedef bspline::Spline<double, K> spline_t;
    std::vector<spline_t> splines;
    // Vector<double> ns; //Normalization constants
    std::vector<size_t> sp_indices; //Unique spline index

    double rmin,rmax; //This might be needed for creating fit basis sets.
    // size_t l;
    // int m;
    
    IrrepIEClient(size_t Ngrid, double rmin, double rmax, size_t l);
    IrrepIEClient(size_t Ngrid, double rmin, double rmax, size_t l,const std::vector<int>& ml);
    virtual ~IrrepIEClient();
    
    virtual size_t size() const {return splines.size();}
    auto indices() const {return ns.indices();}
    auto indices(size_t i) const {return ns.indices(i);}
    const spline_t& operator()(index_t i) const {assert(i>0);return splines[i-1];}

    typedef std::tuple<int,const spline_t*,double> bf_tuple;
    bf_tuple tuple(index_t i) const {return std::make_tuple(l,&splines[i-1],ns(i));}
    GLCache* itsGL;

private:
    IrrepIEClient(const IrrepIEClient<K>&) =delete;
    IrrepIEClient<K>& operator=(const IrrepIEClient<K>&) =delete;
    std::vector<double> MakeLogKnots(size_t NGrid, double rmin, double rmax);
};
} //namespace

