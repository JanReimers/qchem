// File: PolarizedCDImp.H  Implementation for the polarize charge density.
#ifndef _PolarizedCDImp_H_
#define _PolarizedCDImp_H_



#include <ChargeDensity/ChargeDensity.H>

//---------------------------------------------------------------------------------------
//
//  Store spin and spin down a ChargeDensity*'s to allow polymorphism.
//  All member functions just return the unpolarized answer.
//
class Polarized_CDImp
    : public virtual Polarized_CD
{
public:

    Polarized_CDImp(); // No UT coverage
    Polarized_CDImp(DM_CD* up,DM_CD* down);
    ~Polarized_CDImp();

          DM_CD* GetChargeDensity(const Spin&)      ;
    const DM_CD* GetChargeDensity(const Spin&) const;
    
private:
    DM_CD* itsSpinUpCD;
    DM_CD* itsSpinDownCD;
};




#endif //_PolarizedCDImp_H_
