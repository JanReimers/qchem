// File: ElectronContainer.H  Abstract interface for electron containers.
#ifndef _ElectronContainer_H_
#define _ElectronContainer_H_



#include "oml/imp/stream.h"

class ElectronContainer
    : public virtual StreamableObject
{
public:
    virtual ~ElectronContainer() {};
    virtual bool   IsOccupied   (       ) const=0;
    virtual double GetOccupation(       ) const=0;
    virtual double GetSpin      (       ) const=0;
    virtual void   Empty        (       )      =0;
    virtual void   SetOccupation(double )      =0;
    virtual int    GetDegeneracy(       ) const=0;
};

#endif //_ElectronContainer_H_
