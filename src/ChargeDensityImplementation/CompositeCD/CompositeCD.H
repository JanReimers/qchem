// File: CompositeCD.H  Composite charged density, no fitting.
#ifndef _CompositeCD_H_
#define _CompositeCD_H_



#include "ChargeDensity.H"
#include "FunctionsImp/DumpScalarData.H"
#include "Misc/ptr_vector1.h"

//--------------------------------------------------------------------------
//
//  Full charge density represented Compositely as sum of density matricies.
//
class CompositeCD
    : public virtual ChargeDensity
    , private DumpScalarData
{
public:
    CompositeCD();
    void Insert(ChargeDensity*);

    virtual SMat   GetOverlap  (const BasisSet*) const;
    virtual SMat   GetRepulsion(const BasisSet*) const;
    virtual SMat   GetExchange (const BasisSet*) const;

    virtual double GetEnergy(const HamiltonianTerm*) const;

    virtual double GetTotalCharge      (                     ) const;

    virtual void   InjectOverlaps  (FittedFunction*, const BasisSet* theFitBasisSet) const;
    virtual void   InjectRepulsions(FittedFunction*, const BasisSet* theFitBasisSet) const;

    virtual void   ReScale      (double factor              )      ;  //Ro *= factor
    virtual void   ShiftOrigin  (const RVec3&               )      ;
    virtual void   MixIn        (const ChargeDensity&,double)      ;  //this = (1-c)*this + c*that.
    virtual double GetChangeFrom(const ChargeDensity&       ) const;  //MaxAbs(delta density matrix)

    virtual double operator()(const RVec3&) const;
    virtual Vec3   Gradient  (const RVec3&) const;

    virtual std::ostream&       Write(std::ostream&) const;
    virtual std::istream&       Read (std::istream&)      ;

private:
    CompositeCD(const CompositeCD&);
    virtual void  Eval(const Mesh&, Vec&) const;

    optr_vector1<ChargeDensity*> itsCDs;
};

#endif //_CompositeCD_H_
