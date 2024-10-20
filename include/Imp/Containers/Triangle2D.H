// File: Triangle.H  A triangular data structure.
#ifndef _Triangle_H_
#define _Triangle_H_

#include <vector>

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

#endif //_Triangle_H_
