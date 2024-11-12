// File: EnergyLevel.H  Energy level with degeneracy and orbital list.
#ifndef _EnergyLevel_H_
#define _EnergyLevel_H_

#include <types.H>
#include <Spin.H>

struct EnergyLevel1
{
    EnergyLevel1(double _e, double _occ, double _g,const QuantumNumber& _qn,const Spin& _s, Orbital* _o)
        : e(_e), occ(_occ), degen(_g), qn(_qn), s(_s), orbital(_o) 
        {};

    double               e,occ;
    int                  degen;
    const QuantumNumber& qn;
    Spin                 s;
    Orbital*             orbital;
};


#endif //_EnergyLevel_H_
