// File: Symmetry/YlmImp.C  Magnetic spherical harmonic Y_lm(theta,phi) symmetry
module;
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <cassert>
#include <vector>

module qchem.Symmetry.Ylm;
import qchem.Common.Strings;

using std::cout;
using std::endl;


Ylm_Sym::Ylm_Sym(size_t l, const std::vector<int>& _ml) 
: Yl_Sym(l),  ml(_ml) 
{
    assert(ml.size()>0);
};

size_t Ylm_Sym::SequenceIndex() const //Used for op<
 {
    static size_t start=LMAX+1;  //Start after all the Yl Sequence Indexes, LMAX efined in Yl_Sym
    int mmax=*std::max_element(ml.begin(),ml.end());
    return start+mmax+itsL*(itsL+1);
 }


size_t Ylm_Sym::GetDegeneracy() const
{
    return ml.size(); 
}


inline size_t width(int m) {return m<0 ? 2 : 1;}
std::ostream& Ylm_Sym::Write(std::ostream& os) const
{
    os << SPDFG[itsL] << " ";
    if (ml.size()<2*(size_t)itsL+1)
    {
        os << "[";
        for (auto im:ml)
        {
            size_t w=width(im);
            os << std::setw(w) << im << " ";
        }
        os << "]";

    }
   
    return os;
}
