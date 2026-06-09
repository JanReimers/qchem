// File: Symmetry/Internal/Imp/Yl.C  Non magnetic (m-degenerate) spherical harmonic Y_l(theta,phi) symmetry
module;
#include <iostream>
#include <iomanip>
#include <cassert>
#include <blaze/Math.h>

module qchem.Symmetry.Internal.Spherical;
import qchem.Common.Strings;

using std::cout;
using std::endl;

namespace Symmetry::Internal::Spherical
{

Ylm::Ylm(size_t l, const ivec_t& _mls) 
: Yl(l),  mls(_mls) 
{
    assert(mls.size()>0);
};

size_t Ylm::SequenceIndex() const //Used for op<
 {
    static size_t start=LMAX+1;  //Start after all the Yl Sequence Indexes, LMAX efined in Yl
    // int mmax=*std::max_element(mls.begin(),mls.end()); stl veraion
    int mmax=blaze::max(mls);
    return start+mmax+itsL*(itsL+1);
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

} //namespace