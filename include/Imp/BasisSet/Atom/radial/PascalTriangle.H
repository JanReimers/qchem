// File: Imp/Integrals/PascalTriangle.H

#include "Molecule/PolarizedGaussian/Triangle2D.H"

class PascalTriangle : public Triangle2D
{
public:
    static PascalTriangle thePascalTriangle;
private:
    PascalTriangle(int N);

};
