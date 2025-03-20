// File: IEClient.H
#ifndef _IEClient_H_
#define _IEClient_H_

#include <UniqueID.H> 
#include "oml/imp/index_t.h" //for size_t
//
//  Integral DB and engines only sees this
//
class IrrepIEClient //Client for and Irrep basis set.
: public virtual UniqueID
{
public:
    virtual ~IrrepIEClient() {};
    virtual size_t size() const=0;
};

#endif //_IEClient_H_
