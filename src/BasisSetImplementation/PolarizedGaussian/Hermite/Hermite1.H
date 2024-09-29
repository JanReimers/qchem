// File: Hermite1.H  Class for managing 1 function Hermite coefficients
#ifndef _Hermite1_H_
#define _Hermite1_H_

#include <iosfwd>

#define LMAX 3

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

class Polarization;
typedef double Array2D[LMAX+1][LMAX+1];


class Hermite1
{
public:
    Hermite1();
    Hermite1(double AlphaP, int L);
    virtual ~Hermite1() {};

    double  operator()(const Polarization& P,const Polarization& p) const;

    void    Add  (const Hermite1&, double Scale);
    void    Clear(                               );

    virtual std::ostream&  Write(std::ostream&) const;
    virtual std::istream&  Read (std::istream&)      ;
    virtual Hermite1* Clone(        ) const;

private:
    double Getdef(int N,int n) const;
    int     itsL;
    Array2D def;
};

#endif //Hermite1.H
