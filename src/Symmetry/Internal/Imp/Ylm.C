// File: Symmetry/Internal/Imp/Yl.C  Non magnetic (m-degenerate) spherical harmonic Y_l(theta,phi) symmetry
module;
#include <iostream>
#include <iomanip>
#include <cassert>
#include <algorithm>
#include <blaze/Math.h>

module qchem.Symmetry.Internal.Spherical;
import qchem.Common.Strings;

using std::cout;
using std::endl;

namespace Symmetry::Internal::Spherical
{

Ylm::Ylm(size_t l, const ivec_t& _mls) 
: itsL(l), mls(_mls) 
{
    assert(mls.size()>0);
};

size_t Ylm::SequenceIndex() const //Used for op<
 {
    static size_t start=LMax+1;  //Start after all the Yl Sequence Indexes
    
    // Sort mls to get canonical ordering (lexicographic)
    auto sorted_mls = mls;
    std::sort(sorted_mls.begin(), sorted_mls.end());
    
    // Large offset per l value to ensure no collisions between different l values
    size_t l_offset = start + (size_t)itsL * 100000000UL;
    
    // Encode sorted_mls as a positional number with base (2*l+2)
    // This ensures each unique sorted_mls maps to a unique index
    size_t mls_code = 0;
    for (int ml : sorted_mls) {
        int ml_encoded = ml + (int)itsL + 1;  // Map [-l, +l] to [0, 2*l+1]
        mls_code = mls_code * (2*(int)itsL + 2) + ml_encoded;
    }
    
    return l_offset + mls_code;
 }


size_t Ylm::GetDegeneracy() const
{
    return mls.size(); 
}


inline size_t width(int ml) {return ml<0 ? 2 : 1;}
std::ostream& Ylm::Write(std::ostream& os) const
{
    os << SPDFG[itsL] << " ";
    if (mls.size()<2*(size_t)itsL+1)
    {
        os << "{";
        for (auto ml:mls)
        {
            size_t w=width(ml);
            os << std::setw(w) << ml << " ";
        }
        os << "}";

    }
   
    return os;
}

} //namespace