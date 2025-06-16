// File: PascalTriangle.H

#include "Triangle2D.H"

class PascalTriangle : public Triangle2D
{
public:
    static PascalTriangle thePascalTriangle;
private:
    PascalTriangle(int N);

};
