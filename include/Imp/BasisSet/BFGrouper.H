// File: BFGrouper.H  Group Slater or Gaussian basis functions by unique exponents.
#ifndef _BFGrouper_H_
#define _BFGrouper_H_

#include <vector>
#include <map>

// 
// Keep a list of unique exponents for Group Slater or Gaussian basis functions.
// For each unique exponent also store am index and the maximum l angular momentum used
// for that exponent.  The idea is to share exponents between different l-irreps 
// and work out the radial integral tables up to LMax for each exponent combination. 
// THis class should be working together with the charge distribution caching mechanism.
//  
class BFGrouper
{
protected:
    void Append(double e, size_t l);
    size_t LMax(size_t ia, size_t ib, size_t ic, size_t id) const;
    //! Linear array of unique exponents.
    std::vector<double> unique_esv; 
    //! zero based index for each unique exponent.
    std::vector<size_t> es_indices; 

private:
    std::map<double,size_t> unique_es; //Unique exponents.
    //! For each unique exponent, store the maximum l value.
    std::vector<double> maxls;
};



#endif // _BFGrouper_H_
