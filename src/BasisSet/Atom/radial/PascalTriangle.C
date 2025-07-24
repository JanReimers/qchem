// File: BasisSet/Atom/radial/PascalTriangle.C
module;
#include <vector>
#include <iostream>
#include <cassert>
export module qchem.BasisSet.Atom.Internal.radial.PascalTriangle;

#if DEBUG
#define CHECK(i,j) Check(i,j)
#else
#define CHECK(i,j)
#endif
//
//  Implements an triangle data structure, where N+L+M is always < MaxSum.
//
class Triangle2D
{
public:
    Triangle2D(   );
    Triangle2D(size_t);

    double operator()(size_t row,size_t j) const
    {
        CHECK(row,j);
        return itsData[row*(row+1)/2+j];
    }
    double& operator()(size_t row, size_t j)
    {
        CHECK(row,j);
        return itsData[row*(row+1)/2+j];
    }
 
private:
    void    Check(size_t,size_t) const;

    size_t     N;
    std::vector<double> itsData;
};

#undef CHECK

export class PascalTriangle : public Triangle2D
{
public:
    static PascalTriangle thePascalTriangle;
private:
    PascalTriangle(int N);

};

PascalTriangle PascalTriangle::thePascalTriangle(9);


PascalTriangle::PascalTriangle(int N)
: Triangle2D(N+1)
{
    std::cout << "Initializing Pascal Triangle N=" << N << std::endl;
    (*this)(0,0)=1.0;
    for (int row=1;row<=N;row++)
    {
        (*this)(row,0)=(*this)(row,row)=1.0;
        for (int j=1;j<row;j++)
            (*this)(row,j)=(*this)(row-1,j-1)+(*this)(row-1,j);
    }
}

Triangle2D::Triangle2D()
    : N(-1)
{};

Triangle2D::Triangle2D(size_t theMaxSum)
    : N(theMaxSum)
    , itsData  ((N+1)*(N+2)/2)
{};


void Triangle2D::Check(size_t row,size_t j) const
{
    if(row*(row+1)/2+j >= itsData.size())
    {
        std::cerr << "Indices (" << row << "," << j 
                  << ") exceed the Triangle data structure limits"
                  << std::endl << "N = " << N << std::endl;
        assert(false);
        exit(-1);
    }
    if(row<0 || row>=N || j<0 || j>=N)
    {
        std::cerr << "Illegal indices (" << row << "," << j 
                  << ") are not allowed in Triangle data structures"
                  << std::endl << "N = " << N << std::endl;;
        assert(false);
        exit(-1);
    }   
    if(j>row)
    {
        std::cerr << "j >row indices (" << row << "," << j 
                  << ") are not allowed in Triangle data structures"
                  << std::endl << "N = " << N << std::endl;;
        assert(false);
        exit(-1);
    }
}
