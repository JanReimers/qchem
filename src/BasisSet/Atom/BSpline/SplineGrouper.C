// File: BasisSet/Atom/radial/BSpline/SplineGrouper.C Group spline basis functions by unique rmin values.
module;
#include <bspline/Core.h>
#include <vector>
#include <map>
export module qchem.BasisSet.Atom.Internal.SplineGrouper;
export import qchem.BasisSet.Atom.Internal.ExponentGrouper;
import qchem.Basisset.Atom.radial.BSpline.GLQuadrature;

// 
// We can use this class for splines by using the support window {rmin,rmax}.
//
// Keep a list of unique splines for a group BSpline irrep basis functions.
// For each unique rmin,rmax also store an index and the maximum l angular momentum used
// for that spline.  The idea is to share splines between different irreps 
// and work out the radial integral tables up to LMax for each 4-way spline combination. 
// This class should be working together with the charge distribution caching mechanism.
//  
template <class T, size_t K> struct cmpSplines1 {
    bool operator()(const bspline::Spline<T,K>& a, const bspline::Spline<T,K>& b) const
    {
        double af=a.getSupport().front(), bf=b.getSupport().front();
        if (af!=bf) return af<bf;
        return a.getSupport().back()<b.getSupport().back();
    }
};

export template <size_t K> class SplineGrouper : public Grouper
{
    typedef bspline::Spline<double, K> spline_t;
public:
    //! Returns the unique (across all Irrep basis sets) index for this rmin.
    size_t Insert(const spline_t&,size_t l); 
    //! Linear array of unique exponents.
    std::vector<spline_t> unique_spv; 
    std::map<size_t,const GLCache*> itsGLs;
private: 
    //Use this map instead of Grouper::unique_es.
    std::map<spline_t,size_t,cmpSplines1<double,K>> unique_sp; //Unique splines.
};

