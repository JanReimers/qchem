// File: BasisSet/Atom/radial/ExponentGrouper.C Group Slater or Gaussian basis functions by unique exponents.
module;
#include <vector>
#include <map>
export module qchem.BasisSet.Atom.Internal.ExponentGrouper;
// 
// Keep a list of unique exponents for a group Slater or Gaussian irrep basis functions.
// For each unique exponent also store an index and the maximum l angular momentum used
// for that exponent.  The idea is to share exponents between different irreps 
// and work out the radial integral tables up to LMax for each exponent combination. 
// THis class should be working together with the charge distribution caching mechanism.
//
// We can also use this class for splines by using the support window {rmin,rmax}.
// For now just assume if two spline have the same rmin then they also have the same rmax
// So just pass in rmin instead of the exponent.
//  
export class Grouper
{
public:
    size_t LMax(size_t ia, size_t ib, size_t ic, size_t id) const;
protected:
    //! For each unique exponent, store the maximum l value.
    std::vector<size_t> maxls; 
};

export class ExponentGrouper : public Grouper
{
public:
    //! Returns the unique (across all Irrep basis sets) index for this exponent.
    size_t Insert(double exponent,size_t l); 
    //! Linear array of unique exponents.
    std::vector<double> unique_esv; 
private:
    std::map<double,size_t> unique_es; //Unique exponents or rmins for splines.
};

