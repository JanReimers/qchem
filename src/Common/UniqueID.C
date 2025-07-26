// File: UniqueID.C  Anything derived from this will have a unique ID.
export module Common.UniqueID;

export class UniqueID 
{
public:
    virtual ~UniqueID() {};
    typedef int IDtype;
    virtual IDtype GetID() const=0;
};
