// File: UniqueID.H  Anything derived from this will have a unique ID.
#ifndef _UniqueID_H_
#define _UniqueID_H_

class UniqueID 
{
public:
    virtual ~UniqueID() {};
    typedef int IDtype;
    virtual IDtype GetID() const=0;
};

#endif //_UniqueID_H_