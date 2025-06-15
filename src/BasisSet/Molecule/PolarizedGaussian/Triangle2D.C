// File: Triangle.C  A triangular data structure.
#include "Imp/BasisSet/Molecule/PolarizedGaussian/Triangle2D.H"
#include "Common/stl_io.h"
#include <cassert>

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

