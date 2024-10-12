// File: ExactIrrepCD.H  Exact charged density for ONE irreducable representation basis set.
#ifndef _ExactIrrepCD_H_
#define _ExactIrrepCD_H_



#include "BasisSet.H"
#include "ChargeDensity.H"
#include "Misc/Spin.H"
#include "oml/smatrix.h"

template <class T> class TIrrepBasisSetCommon;

//------------------------------------------------------------------------------------
//
//  This maintains the exact charge density represented by the density matrix
//  of one irreducable representation.  The full charge density will in general
//  be a summation of these guys.
//
template <class T> class ExactIrrepCD
    : public virtual ChargeDensity
{
public:
    typedef Matrix<T>  DenMat;
    typedef SMatrix<T> DenSMat; //Type for the density matrix.

    ExactIrrepCD();
    ExactIrrepCD(const DenSMat&,const IrrepBasisSet*, const Spin&);

    virtual SMat   GetOverlap  (const IrrepBasisSet*) const;
    virtual SMat   GetRepulsion(const IrrepBasisSet*) const;
    virtual SMat   GetExchange (const IrrepBasisSet*) const;

    virtual double GetEnergy     (const HamiltonianTerm*) const;
    virtual double GetTotalCharge(                ) const;

    virtual Vector<double> GetRepulsions(const IrrepBasisSet* theFitBasisSet) const;

    virtual void   ReScale      (double factor              )      ;
    virtual void   ShiftOrigin  (const RVec3&               )      ;
    virtual void   MixIn        (const ChargeDensity&,double)      ;  //this = (1-c)*this + c*that.
    virtual double GetChangeFrom(const ChargeDensity&       ) const;  //MaxAbs(delta density matrix)

    virtual double operator()(const RVec3&) const;
    virtual Vec3   Gradient  (const RVec3&) const;

    virtual std::ostream&       Write(std::ostream&) const;
    virtual std::istream&       Read (std::istream&)      ;

private:

    virtual void  Eval(const Mesh&, Vec&) const;

    DenSMat                  itsDensityMatrix;
    const IrrepBasisSet*     itsBasisSet;
    const TIrrepBasisSet<T>* itsCastedBasisSet;
    Spin                     itsSpin;
};

#endif //_ExactIrrepCD_H_
