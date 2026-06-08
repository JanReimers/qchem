// File: Symmetry/YlmImp.C  Magnetic spherical harmonic Y_lm(theta,phi) symmetry
module;
// #include <algorithm>
#include <iostream>
#include <iomanip>
#include <cassert>
// #include <vector>
#include <blaze/Math.h>
// #include <blaze/util/algorithms/Max.h>

module qchem.Symmetry.Ylm;
import qchem.Common.Strings;

using std::cout;
using std::endl;


Ylm_Sym::Ylm_Sym(size_t l, const ivec_t& _mls) 
: Yl_Sym(l),  mls(_mls) 
{
    assert(mls.size()>0);
};

size_t Ylm_Sym::SequenceIndex() const //Used for op<
 {
    static size_t start=LMAX+1;  //Start after all the Yl Sequence Indexes, LMAX efined in Yl_Sym
    // int mmax=*std::max_element(mls.begin(),mls.end()); stl veraion
    int mmax=blaze::max(mls);
    return start+mmax+itsL*(itsL+1);
 }


size_t Ylm_Sym::GetDegeneracy() const
{
    return mls.size(); 
}


inline size_t width(int ml) {return ml<0 ? 2 : 1;}
std::ostream& Ylm_Sym::Write(std::ostream& os) const
{
    os << SPDFG[itsL] << " ";
    if (mls.size()<2*(size_t)itsL+1)
    {
        os << "[";
        for (auto ml:mls)
        {
            size_t w=width(ml);
            os << std::setw(w) << ml << " ";
        }
        os << "]";

    }
   
    return os;
}
