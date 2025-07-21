// File: BFGrouper.H  Group Slater or Gaussian basis functions by unique exponents.
module;
#include <vector>
#include <map>
export module qchem.BasisSet.Atom.BFGrouper;
export import qchem.BasisSet.Atom.IEClient;

// 
// Keep a list of unique exponents for Group Slater or Gaussian basis functions.
// For each unique exponent also store am index and the maximum l angular momentum used
// for that exponent.  The idea is to share exponents between different l-irreps 
// and work out the radial integral tables up to LMax for each exponent combination. 
// THis class should be working together with the charge distribution caching mechanism.
//  

export class BFGrouper
{
protected:
    void Append(AtomIrrepIEClient*);
    size_t LMax(size_t ia, size_t ib, size_t ic, size_t id) const;
    //! Linear array of unique exponents.
    std::vector<double> unique_esv; 
private:
    std::map<double,size_t> unique_es; //Unique exponents.
    //! For each unique exponent, store the maximum l value.
    std::vector<double> maxls;
};
