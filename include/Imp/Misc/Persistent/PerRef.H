// File: PerRef.H  Reference to a persistent object.
#ifndef _PerRef_H_
#define _PerRef_H_



#include "Misc/pmstream.h"
#include "Misc/UniqueID/UniqueID.H"

#if DEBUG
#define CHECK Check()
#else
#define CHECK
#endif

template <class T> class PerRef
    : public virtual PMStreamableObject
{
public:
    PerRef(      T*);
    PerRef(        ) : itsPointer(0), itsID(0) {};

    const T*   operator->(               ) const
    {
        CHECK;
        return  itsPointer;
    }
    T*   operator->(               )
    {
        CHECK;
        return  itsPointer;
    }
    const T&   operator* (               ) const
    {
        CHECK;
        return *itsPointer;
    }
    T&   operator* (               )
    {
        CHECK;
        return *itsPointer;
    }
    operator const T*    (               ) const
    {
        CHECK;
        return  itsPointer;
    }
    operator       T*    (               )
    {
        CHECK;
        return  itsPointer;
    }
    const bool operator==(const PerRef& i) const
    {
        CHECK;
        return  itsPointer==i.itsPointer && itsID==i.itsID;
    }

    virtual std::ostream& Write     (std::ostream&) const;
    virtual std::istream& Read      (std::istream&)      ;
    virtual PerRef*  Clone     (        ) const;
    static  PerRef*  Factory   (std::istream&)      ;

private:
    void  Check() const;

    T*               itsPointer;
    UniqueID::IDtype itsID;
};

#undef CHECK

#endif
