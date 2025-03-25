// File: ExactIrrepCD.H  Exact charged density for ONE irreducable representation basis set.
#ifndef _ExactIrrepCD_H_
#define _ExactIrrepCD_H_



#include <ChargeDensity.H>
#include <Spin.H>
#include "oml/smatrix.h"

//------------------------------------------------------------------------------------
//
//  This maintains the exact charge density represented by the density matrix
//  of one irreducable representation.  The full charge density will in general
//  be a summation of these guys.
//
template <class T> class IrrepCD
    : public virtual Exact_CD
{
public:
    typedef Matrix<T>  DenMat;
    typedef SMatrix<T> DenSMat; //Type for the density matrix.
    typedef DensityFFClient::RVec RVec;

    IrrepCD();
    IrrepCD(const DenSMat&,const TOrbital_IBS<T>*, const Spin&);

    virtual SMat   GetRepulsion(const TOrbital_IBS<double>*) const;
    virtual SMat   GetExchange (const TOrbital_IBS<double>*) const;
    virtual RVec   GetRepulsion3C(const Fit_IBS*) const;

    virtual double GetEnergy     (const HamiltonianTerm*) const;
    virtual double GetTotalCharge(                      ) const;

    virtual void   ReScale      (double factor         )      ; // No UT coverage
    virtual void   ShiftOrigin  (const RVec3&          )      ; // No UT coverage
    virtual void   MixIn        (const Exact_CD&,double)      ;  //this = (1-c)*this + c*that.
    virtual double GetChangeFrom(const Exact_CD&       ) const;  //MaxAbs(delta density matrix)

    virtual double operator()(const RVec3&) const;
    virtual Vec3   Gradient  (const RVec3&) const; // No UT coverage

    virtual std::ostream&       Write(std::ostream&) const;
    virtual std::istream&       Read (std::istream&)      ;

private:
    bool IsZero() const;
    SMat ZeroM(size_t N) const;
    RVec ZeroV(size_t N) const;

    DenSMat                itsDensityMatrix;
    const TOrbital_IBS<T>* itsBasisSet;
    Spin                   itsSpin;
};

#endif //_ExactIrrepCD_H_
