// File: IEClient.H
#ifndef _IEClient_H_
#define _IEClient_H_

#include "Imp/Misc/UniqueID/UniqueID.H" //Just need ID type.
//
//  Integral DB and engines only sees this
//
class IrrepIEClient //Client for and Irrep basis set.
{
public:
    virtual ~IrrepIEClient() {};
    virtual UniqueID::IDtype GetID() const=0;
    virtual size_t size() const=0;
};

#endif //_IEClient_H_
