// File: BasisSet/Atom/Evaluators/Internal/ExponentGrouper.C Group Slater or Gaussian basis functions by unique exponents.
module;
#include <vector>
#include <map>
#include "forward.H"
export module qchem.BasisSet.Atom.Evaluators.Internal.Grouper;

namespace qchem {
// 
// Keep a list of unique exponents for a group Slater or Gaussian irrep basis functions.
// For each unique exponent also store an index and the maximum l angular momentum used
// for that exponent.  The idea is to share exponents between different irreps 
// and work out the radial integral tables up to LMax for each exponent combination. 
// THis class should be working together with the charge distribution caching mechanism.
//
// The spline analogue (SplineGrouper) keys on a spline's full identity -- its knot vector plus support
// window -- so two splines from different grids never alias even when they share a single "BSpline<K>"
// Cache4.  (It deliberately does NOT reduce a spline to a scalar like rmin: that was lossy and only worked
// while each grid had its own Cache4 keyed by the whole grid string.)
//
export class Grouper
{
public:
    size_t LMax(size_t ia, size_t ib, size_t ic, size_t id) const;
protected:
    friend class ::Cache4Tests;
    //! For each unique exponent, store the maximum l value.
    std::vector<size_t> maxls; 
};

export class ExponentGrouper : public Grouper
{
public:
    ~ExponentGrouper() {}; //g++ 15.2 BUG compiler generated version will not instance std::map::~map()
    //! Returns the unique (across all Irrep basis sets) index for this exponent.
    size_t Insert(double exponent,size_t l); 
    //! Linear array of unique exponents.
    std::vector<double> unique_esv; 
private:
    friend class ::Cache4Tests;
    std::map<double,size_t> unique_es; //Unique exponents or rmins for splines.
};


} // namespace qchem