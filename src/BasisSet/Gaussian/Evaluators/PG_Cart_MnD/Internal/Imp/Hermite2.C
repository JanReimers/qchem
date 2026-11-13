// File: BasisSet/Gaussian/Evaluators/PG_Cart_MnD/Internal/Imp/Hermite2.C  Class for managing 2 function Hermite coefficients
module;

#include <iostream>
#include <cassert>
#include <iomanip>
#include <vector>

module qchem.BasisSet.Gaussian.Evaluators.PG_Cart_MnD.Internal.Hermite;
import qchem.Reporting;   // FormatGuard
import qchem.BasisSet.Gaussian.Evaluators.PG_Cart_MnD.Polarization;

std::ostream& operator << (std::ostream& os, const std::vector<double>& v)
{
    for (auto i:v) os << i << " ";
    return os;
}

namespace qchem::BasisSet::Gaussian::Evaluators::PG_Cart_MnD
{
//----------------------------------------------------------------------------------------
//
//  None of this will make any sense unless you read McMurchie and Davidson,
//  J. Comp. Phys. 26, 218 (1978).
//

/*
 The storage here is non-trivial:
 I will label the elements as d(N,n1,n2) with 0 <= n1+n2 = n < LMax, and 0 < N <= n.
 d(0,n1,n2) looks like this up to n=2:

    1         B       B**2+1/2a .
    A         AB+1/2a .         .
    A**2+1/2a .       .         .
    .         .       .         .


 d(1,n1,n2) looks like

    0     1/2a     2B/2a .
    1/2a  (A+B)/2a .     .
    2A/2a .        .     .
    .     .        .     .

 d(2,n1,n2) looks like

    0         0         1/(2a)**2 .
    0         1/(2a)**2 .         .
    1/(2a)**2 .         .         .
    .         .         .         .

 You can see the all the 0 and 1/(2a)**N elements need not be stored.   And all the zeros
also need not be stored.  The elements are stored using the fancy lookup table below to transform
(N,na,nb) ---> array index.  ANother lookup table is used to get the required size for a given
LA,LB pair.
   The N arrays are filled in the following order:  As you move from top left to bottom right
 n=n1+n2 increases.  Start by loading the n=1 diagonal for n2=0 to 1, in the N=0 array.  N>0
 does not need this row filled.  The move on to n=2, and do n2=0 to 2 for N=0 and 1.  Then n=3
 for N=1 to n-1, and n2=0 to n+1.
*/

using std::cout;
using std::endl;
//----------------------------------------------------------------------------------------




inline bool near(double a, double b) 
{
    static double eps=1e-13;
    double err=fabs(a)<1 ? fabs(a-b) : fabs(a-b)/fabs(a) ;
    if (err>=eps)
    {
        qchem::report::FormatGuard keep(cout);
        cout << std::setprecision(14) << a << " " << b << " " << err << endl;
    }
    return err<eps;
}

inline bool near(const rvec3_t& a, const rvec3_t& b) 
{
    return near(a.x,b.x) && near(a.y,b.y) &&  near(a.z,b.z);
}

// The constructor, destructor, Assign and operator() are INLINE in the module interface (Enzyme differentiates through them).
std::ostream& Hermite2::Write(std::ostream& os) const
{
    qchem::report::FormatGuard keep(os);
    os << ": LA=" << LA << ", LB = " << LB << std::endl;
    os.precision(4);
    os.width(8);
    os.setf(std::ios::fixed,std::ios::floatfield);
    for (int N=0; N<=LA+LB; N++)
    {
        os << "---------------------------------------------------------------------" << std::endl;
        os << "N = " << N << std::endl;
        for (int na=0; na<=LA; na++)
        {
            for (int nb=0; nb<=LB; nb++)  os << d[GetIndex(N,na,nb)] << " ";
            os << std::endl;
        }
    }
    for (int N=0; N<=LA+LB; N++)
    {
        os << "---------------------------------------------------------------------" << std::endl;
        os << "L = " << N << std::endl;
        for (int na=0; na<=LA; na++)
        {
            for (int nb=0; nb<=LB; nb++)  os << e[GetIndex(N,na,nb)] << " ";
            os << std::endl;
        }
    }
    for (int N=0; N<=LA+LB; N++)
    {
        os << "---------------------------------------------------------------------" << std::endl;
        os << "M = " << N << std::endl;
        for (int na=0; na<=LA; na++)
        {
            for (int nb=0; nb<=LB; nb++)  os << f[GetIndex(N,na,nb)] << " ";
            os << std::endl;
        }
    }
    return os;
}


Hermite2* Hermite2::Clone() const
{
    return new Hermite2(*this);
}


} //namespace qchem::BasisSet::Gaussian::Evaluators::PG_Cart_MnD
