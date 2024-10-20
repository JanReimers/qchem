// File: UniqueID.H  Anything derived from this will have a unique ID.
#ifndef _UniqueID_H_
#define _UniqueID_H_



#include "Imp/Misc/pmstream.h"

//#############################################################################
//
//  This class allows objects to be pickled with unique ID's.  The problem
//  is in assigning unused ID's to new objects.  Within a single heavy weight
//  process a statice counter is used.  When ID'd objects are destroyed
//  the current count (NectID) is saved in a file on disk.  When a new
//  process is started this file will be used to initialize the NextID.
//
class UniqueID :
    public virtual PMStreamableObject
{
public:
    typedef int IDtype;

    UniqueID();
    UniqueID(const UniqueID&);
    virtual ~UniqueID();
    UniqueID& operator=(const UniqueID&);

    virtual std::ostream& Write(std::ostream&) const=0;
    virtual std::istream& Read (std::istream&)      =0;

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

#endif
