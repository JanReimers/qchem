// FIle: SCFAccelerator_Null.H  A simple pass through accerlator proxy that does no acceleration.
#ifndef _SCFAccelerator_Null_H_
#define _SCFAccelerator_Null_H_

#include "SCFAccelerator.H"

#include "oml/smatrix.h"
class SCFIrrepAccelerator__Null : public virtual SCFIrrepAccelerator
{
public:
    SCFIrrepAccelerator__Null(const Irrep_QNs&) {};
    virtual ~SCFIrrepAccelerator__Null() {};
    virtual void Init(const LASolver<double>*); 
    virtual void UseFD(const SMat& F, const SMat& DPrime);
    virtual SMat Project();
private:
    const LASolver<double>*   itsLaSolver; //Knows the ortho transform
    SMat   itsFPrime;
};



#endif //_SCFAccelerator_Null_H_