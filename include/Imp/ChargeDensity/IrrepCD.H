// File: ExactIrrepCD.H  Exact charged density for ONE irreducable representation basis set.
#ifndef _ExactIrrepCD_H_
#define _ExactIrrepCD_H_



#include <BasisSet.H>
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
    : public virtual ChargeDensity
{
public:
    typedef Matrix<T>  DenMat;
    typedef SMatrix<T> DenSMat; //Type for the density matrix.
    typedef DensityFFClient::RVec RVec;

    IrrepCD();
    IrrepCD(const DenSMat&,const Orbital_IBS<T>*, const Spin&);

    virtual SMat   GetRepulsion(const Orbital_IBS<double>*) const;
    virtual SMat   GetExchange (const Orbital_IBS<double>*) const;
    virtual RVec   GetRepulsion3C(const IrrepBasisSet* theFitBasisSet) const;

    virtual double GetEnergy     (const HamiltonianTerm*) const;
    virtual double GetTotalCharge(                ) const;

    virtual void   ReScale      (double factor              )      ;
    virtual void   ShiftOrigin  (const RVec3&               )      ;
    virtual void   MixIn        (const ChargeDensity&,double)      ;  //this = (1-c)*this + c*that.
    virtual double GetChangeFrom(const ChargeDensity&       ) const;  //MaxAbs(delta density matrix)

    virtual double operator()(const RVec3&) const;
    virtual Vec3   Gradient  (const RVec3&) const;

    virtual std::ostream&       Write(std::ostream&) const;
    virtual std::istream&       Read (std::istream&)      ;

private:
    bool IsZero() const;
    SMat ZeroM(const IrrepBasisSet* bs_ab) const;
    RVec ZeroV(const IrrepBasisSet* bs_ab) const;

    DenSMat                  itsDensityMatrix;
    const Orbital_IBS<T>* itsBasisSet;
    Spin                     itsSpin;
};

#endif //_ExactIrrepCD_H_
