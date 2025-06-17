// File: ExactIrrepCD.H  Exact charged density for ONE irreducable representation basis set.
#ifndef _ExactIrrepCD_H_
#define _ExactIrrepCD_H_



#include <ChargeDensity/ChargeDensity.H>
#include <Symmetry/Spin.H>
#include <Symmetry/Irrep_QNs.H>
#include "oml/smatrix.h"

//------------------------------------------------------------------------------------
//
//  This maintains the exact charge density represented by the density matrix
//  of one irreducable representation.  The full charge density will in general
//  be a summation of these guys.
//
template <class T> class IrrepCD
    : public virtual DM_CD
{
public:
    typedef Matrix<T>  DenMat;
    typedef SMatrix<T> DenSMat; //Type for the density matrix.
    typedef DensityFFClient::RVec RVec;

    IrrepCD();
    IrrepCD(const DenSMat&,const TOrbital_IBS<T>*, Irrep_QNs);

    virtual SMat   GetRepulsion(const TOrbital_HF_IBS<double>*) const;
    virtual SMat   GetExchange (const TOrbital_HF_IBS<double>*) const;
    virtual RVec   GetRepulsion3C(const Fit_IBS*) const;

    virtual double DM_Contract(const Static_CC*) const;
    virtual double DM_Contract(const Dynamic_CC*,const DM_CD*) const;
    virtual double GetTotalCharge(                      ) const;

    virtual void   ReScale      (double factor         )      ; // No UT coverage
    virtual void   ShiftOrigin  (const RVec3&          )      ; // No UT coverage
    virtual void   MixIn        (const DM_CD&,double)      ;  //this = (1-c)*this + c*that.
    virtual double GetChangeFrom(const DM_CD&       ) const;  //MaxAbs(delta density matrix)

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
    Irrep_QNs              itsIrrep;
};

#endif //_ExactIrrepCD_H_
