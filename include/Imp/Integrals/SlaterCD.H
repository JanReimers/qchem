// File: SlaterCD.H  4 electron Charge distribution of Slater orbitals. 
#ifndef _SlaterCD_H_
#define _SlaterCD_H_

#include "Imp/Integrals/Cache4.H"
#include "oml/matrix.h"
//
// Store exponents and derivative tables for a 4 electron charge distribution ab|cd
// formed from a product of slater orbitals.
// eab = ea + eb
// ecd = ec + ed
//
class SlaterCD : public virtual Cacheable
{
public:
    SlaterCD(double eab, double ecd, size_t LMax);
    double Coulomb_R0(size_t la,size_t lc) const; //R_0(la,la,lc,lc);
    Vector<double> Coulomb_Rk(size_t la,size_t lc) const; //R_k(la,la,lc,lc);
    Vector<double> ExchangeRk(size_t la,size_t lb) const; //R_k(la,lb,la,lb);
    
private:
    static double fk(double a, double ab, size_t k,size_t n);
    
    double eab, ecd;
    size_t LMax;
    Matrix<double> Iab,Icd; //Derivative tables.
};

#endif // _SlaterCD_H_
