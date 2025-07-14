// File: UniqueIDImp.H  Anything derived from this will have a unique ID.
#ifndef _UniqueIDImp_H_
#define _UniqueIDImp_H_


#include <Common/UniqueID.H>
#include "Common/pmstream.h"
#include <iostream>
//#############################################################################
//
//  This class allows objects to be pickled with unique ID's.  The problem
//  is in assigning unused ID's to new objects.  Within a single heavy weight
//  process a statice counter is used.  When ID'd objects are destroyed
//  the current count (NectID) is saved in a file on disk.  When a new
//  process is started this file will be used to initialize the NextID.
//
class UniqueIDImp 
    : public virtual UniqueID
    , public virtual PMStreamableObject
{
public:
    typedef int IDtype;

    UniqueIDImp();
    UniqueIDImp(const UniqueIDImp&);
    ~UniqueIDImp();
    
    UniqueID& operator=(const UniqueID&);

    virtual std::ostream& Write(std::ostream&) const;

    IDtype GetID() const
    {
        return itsID;
    }
    IDtype NewID()
    {
        return itsID=GetNextID();
    }
    static
    IDtype GetMaxID()
    {
        return MaxID;
    }

private:
    static IDtype GetNextID();
    IDtype itsID;
    static       IDtype NextID;
    static const IDtype MaxID;
};

#endif //_UniqueIDImp_H_
