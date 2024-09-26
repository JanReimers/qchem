// File: IDRef.H  Reference to an ID'd object.
#ifndef _IDRef_H_
#define _IDRef_H_



#include "Misc/pmstream.h"
#include "Misc/UniqueID/UniqueID.H"

#if DEBUG
#define CHECK Check()
#else
#define CHECK
#endif

template <class T> class IDRef
    : public virtual PMStreamableObject
{
public:
    IDRef(T*);
    IDRef(  ) : itsPointer(0), itsID(0) {};
    void FixUpPointer(T*);

    const T*   operator->() const
    {
        CHECK;
        return  itsPointer;
    }
    T*   operator->()
    {
        CHECK;
        return  itsPointer;
    }
    const T&   operator* () const
    {
        CHECK;
        return *itsPointer;
    }
    T&   operator* ()
    {
        CHECK;
        return *itsPointer;
    }
    operator const T*    () const
    {
        CHECK;
        return  itsPointer;
    }
    operator       T*    ()
    {
        CHECK;
        return  itsPointer;
    }

    const bool operator==(const IDRef& i) const
    {
        return itsID==i.itsID && itsID!=0;
    }

    virtual std::ostream& Write     (std::ostream&) const;
    virtual std::istream& Read      (std::istream&)      ;
    virtual IDRef*   Clone     (        ) const;
    static  IDRef*   Factory   (std::istream&)      ;

private:
//    typedef PMStreamableObject Base;
    void  Check() const;

    T*               itsPointer;
    UniqueID::IDtype itsID;
};

#undef CHECK

#endif
