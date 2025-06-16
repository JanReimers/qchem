// File: Integrals/PascalTriangle.C

#include "Atom/radial/PascalTriangle.H"
#include <iostream>

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
