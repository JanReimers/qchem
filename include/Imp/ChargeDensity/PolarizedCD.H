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
class Polarized_Exact_CDImp
    : public virtual Polarized_Exact_CD
{
public:

    Polarized_Exact_CDImp(); // No UT coverage
    Polarized_Exact_CDImp(Exact_CD* up,Exact_CD* down);
    ~Polarized_Exact_CDImp();

          Exact_CD* GetChargeDensity(const Spin&)      ;
    const Exact_CD* GetChargeDensity(const Spin&) const;

private:
    Exact_CD* itsSpinUpCD;
    Exact_CD* itsSpinDownCD;
};

#endif //_PolarizedCDImp_H_
