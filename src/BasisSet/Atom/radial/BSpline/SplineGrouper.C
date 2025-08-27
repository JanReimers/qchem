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
// For now just assume if two spline have the same rmin then they also have the same rmax
// So just pass in rmin instead of the exponent.
//
// Keep a list of unique splines for a group BSpline irrep basis functions.
// For each unique rmin also store an index and the maximum l angular momentum used
// for that spline.  The idea is to share splines between different irreps 
// and work out the radial integral tables up to LMax for each 4-way spline combination. 
// This class should be working together with the charge distribution caching mechanism.
//  
export template <size_t K> class SplineGrouper : public Grouper
{
    typedef bspline::Spline<double, K> spline_t;
public:
    //! Returns the unique (across all Irrep basis sets) index for this rmin.
    size_t Insert(const spline_t&,size_t l); 
    //! Linear array of unique exponents.
    std::vector<spline_t> unique_spv; 
    std::map<size_t,const GLCache*> itsGLs;
};

