// File: PolarizedGaussian/MnD/Hermite1.C  Class for managing 1 function Hermite coefficients
module;

#include <iosfwd>

#define LMAX 3
export module qchem.BasisSet.Molecule.PolarizedGaussian.Internal.MnD.Hermite1;
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.Polarization;

//
//  Encapsulate calculation and storage of Hermite expansion coefficients
//
//   n  l  m
//  d  e  f
//   N  L  M
//
//  The op() function returns the product of all three.  This is more efficient
//  because this function can do fast checks to see if any of the three are
//  zero before doing a more expensive lookup.
//

typedef double Array2D[LMAX+1][LMAX+1];

export namespace PolarizedGaussian
{

class Hermite1
{
public:
    Hermite1();
    Hermite1(double AlphaP, int L);
 
    double  operator()(const Polarization& P,const Polarization& p) const;

    void    Add  (const Hermite1&, double Scale);
    void    Clear();

    // std::ostream&  Write(std::ostream&) const; //For debugging.

private:
    double Getdef(int N,int n) const;
    int     itsL;
    Array2D def;
};

} //namespace PolarizedGaussian

