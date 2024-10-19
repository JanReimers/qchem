// File: PolarizedCDImp.H  Implementation for the polarize charge density.
#ifndef _PolarizedCDImp_H_
#define _PolarizedCDImp_H_



#include <ChargeDensity.H>
#include <Spin.H>

//---------------------------------------------------------------------------------------
//
//  Store spin and spin down a ChargeDensity*'s to allow polymorphism.
//  All member functions just return the unpolarized answer.
//
class PolarizedCDImp
    : public virtual PolarizedCD
{
public:

    PolarizedCDImp();
    PolarizedCDImp(ChargeDensity* up, ChargeDensity* down);
    ~PolarizedCDImp();

    ChargeDensity*       GetChargeDensity(const Spin&)      ;
    const ChargeDensity* GetChargeDensity(const Spin&) const;

    virtual std::ostream&        Write(std::ostream&) const;
    virtual std::istream&        Read (std::istream&)      ;

private:
    ChargeDensity* itsSpinUpCD;
    ChargeDensity* itsSpinDownCD;
};

#endif //_PolarizedCDImp_H_
