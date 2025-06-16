// File: Triangle.C  A triangular data structure.
#include "PolarizedGaussian/Triangle3D.H"
#include "Common/stl_io.h"
#include <cassert>

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
        std::cerr << "Indecies (" << i << "," << j << "," << k
                  << ") exceed the Triangle data structure limits"
                  << std::endl << "MaxSum = " << N << std::endl;
        assert(false);
        exit(-1);
    }
    if(i<0 || j<0 || k<0)
    {
        std::cerr << "Negative indecies (" << i << "," << j << "," << k
                  << ") are not allowed in Triangle data structures"
                  << std::endl << "MaxSum = " << N << std::endl;;
        assert(false);
        exit(-1);
    }
}

std::ostream&  Triangle3D::Write(std::ostream& os) const
{
    if (!StreamableObject::Pretty())
        os << N << " " << itsData << std::endl;
    else
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
    }

    return os;
}

std::istream&  Triangle3D::Read (std::istream& is)
{
    return is >> N >> itsData;
}

