// File: PolarizedCDImplementation.H  Implementation for the polarize charge density category.
#ifndef _PolarizedCDImplementation_H_
#define _PolarizedCDImplementation_H_



#include "ChargeDensity/PolarizedCD.H"
#include "FunctionsImp/DumpScalarData.H"
#include "Misc/Spin.H"

//---------------------------------------------------------------------------------------
//
//  Store spin and spin down a ChargeDensity*'s to allow polymorphism.
//  All member functions just return the unpolarized answer.
//
class PolarizedCDImplementation
    : public virtual PolarizedCD
    , private DumpScalarData
{
public:

    PolarizedCDImplementation();
    PolarizedCDImplementation(ChargeDensity* up, ChargeDensity* down);
    ~PolarizedCDImplementation();

    ChargeDensity* GetChargeDensity(const Spin&)      ;
    const ChargeDensity* GetChargeDensity(const Spin&) const;

    virtual std::ostream&        Write(std::ostream&) const;
    virtual std::istream&        Read (std::istream&)      ;

private:
    ChargeDensity* itsSpinUpCD;
    ChargeDensity* itsSpinDownCD;
};

#endif //_PolarizedCDImplementation_H_
