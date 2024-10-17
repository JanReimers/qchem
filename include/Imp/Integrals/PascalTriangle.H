// File: Imp/Integrals/PascalTriangle.H

#include "Imp/Containers/Triangle2D.H"

class PascalTriangle : public Triangle2D
{
public:
    static PascalTriangle thePascalTriangle;
private:
    PascalTriangle(int N);

};
