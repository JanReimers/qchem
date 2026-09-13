// File: BasisSet/Gaussian/Evaluators/Internal/MnD/Imp/Triangle3D.C  A 3-index triangular data structure.
module;
#include <cassert>
#include <iostream>
#include <stdexcept>   // R2.5b: an out-of-range triangle index THROWS
#include <string>
#include <vector>
module qchem.BasisSet.Gaussian.Evaluators.Internal.MnD.Triangle3D;
import qchem.stl_io;

namespace qchem::BasisSet::Gaussian::Evaluators::Internal::MnD
{

Triangle3D::Triangle3D()
    : N(-1)
{};

Triangle3D::Triangle3D(int theMaxSum)
    : N(theMaxSum)
    , itsData  ((N+1)*(N+2)*(N+3)/6)
{};

Triangle3D& Triangle3D::operator=(const Triangle3D& other)
{
    if(&other!=this)
    {
        N = other.N;
        itsData   = other.itsData;
    }
    return *this;
}

void Triangle3D::Add(const Triangle3D& theT, double theScale)
{
    assert(N<0 || N==theT.N);
    if (N<0)
    {
        N=theT.N;
        itsData.resize(theT.itsData.size(),0.0);
    }
    std::vector<double>::iterator i(itsData.begin());
    std::vector<double>::const_iterator  b(theT.itsData.begin());
    for (; i!=itsData.end()&&b!=itsData.end(); i++,b++) *i+=*b * theScale;
}


void Triangle3D::Check(int i,int j,int k) const
{
    if(i+j+k > N)
    {
        throw std::out_of_range("Triangle3D: indices ("+std::to_string(i)+","+std::to_string(j)+","
            +std::to_string(k)+") sum to "+std::to_string(i+j+k)+", past this triangle's MaxSum of "
            +std::to_string(N)+".  The caller built the table for one angular reach and indexed it for a "
            "larger one -- a composition error, not a recoverable condition.");
    }
    if(i<0 || j<0 || k<0)
    {
        throw std::out_of_range("Triangle3D: negative indices ("+std::to_string(i)+","+std::to_string(j)
            +","+std::to_string(k)+").  A triangle is indexed by Cartesian POWERS, which are >= 0; a "
            "negative one means a recurrence stepped below its base case without taking the branch.");
    }
}

std::ostream&  Triangle3D::Write(std::ostream& os) const
{
    os << "Triangle N=" << N << "data size=" << itsData.size() << std::endl;
    for (int i=0; i<=N; i++)
    {
        os << "i = " << i << " layer:" << std::endl;
        for (int j=0; j<=N-i; j++)
        {
            for (int k=0; k<=N-j-i; k++) os << (*this)(i,j,k) << " ";
            os << std::endl;
        }
    }
    return os;
}

} //namespace


