// File: BFGrouper.H  Group Slater or Gaussian basis functions by unique exponents.
#ifndef _BFGrouper_H_
#define _BFGrouper_H_

#include <vector>
#include <map>

//
// 
// Group Slater or Gaussian basis functions by unique exponents.  When using this
// class there are two sets of indices:
//   1) flat indices which include duplicates
//   2) unique indices which refer to the list of unique exponents.
// 
//  
class BFGrouper
{
public:
    //! zero based index for each unique exponent.
    std::vector<size_t> es_indices; 
    //! Linear array of unique exponents.
    std::vector<double> unique_esv; 
    //! Get all indices with Ls(i)=l
    const std::vector<size_t>& indices(size_t l) const;

protected:
    void Append(double e, size_t l, size_t flat_index);

private:
    std::map<double,size_t> unique_es; //Unique exponents.
    std::map<size_t,std::vector<size_t> > L_indices; // <L,flat indices>
};



#endif // _BFGrouper_H_
