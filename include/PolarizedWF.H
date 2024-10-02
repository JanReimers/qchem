// File: PolarizedWF.H  Interface for a polarized wave function.
#ifndef _PolarizedWF_H_
#define _PolarizedWF_H_



#include "WaveFunction.H"
#include "PolarizedCD.H"

class PolarizedWF
    : public virtual WaveFunction
{
public:
    virtual WaveFunction* GetWaveFunction (const Spin&)=0;
    static  PolarizedWF*  Factory(std::istream&);
};

#endif //_PolarizedWF_H_
