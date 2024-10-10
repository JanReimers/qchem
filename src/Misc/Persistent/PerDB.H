// File: PerDB.H  Manage a global list of ID object *'s.
#ifndef _PerDB_H_
#define _PerDB_H_



#include "Misc/UniqueID/UniqueID.H"
#include "Imp/Containers/ptr_vector.h"

template <class T> class PerRef;
template <class T> class PerList;

template <class T> class PerDB
    : public virtual PMStreamableObject
{
public:
    PerDB();
    ~PerDB() {thePerDB=0;}
    virtual std::ostream& Write(std::ostream&) const;
    virtual std::istream& Read (std::istream&)      ;

    void Dump(std::ostream&) const;
    void Clear();

private:
    static T*   GetPointerFor(UniqueID::IDtype);
    static bool StaticHas    (UniqueID::IDtype);
    static T*   Insert       (      T*        );

    UniqueID::IDtype Find        (const T&        ) const;
    bool             Has         (UniqueID::IDtype) const;
    void             Delete      (UniqueID::IDtype)      ;

//  virtual PerDB*   Clone(        ) const;

    friend class PerRef <T>;
    friend class PerList<T>;

    optr_vector1<T*>  itsPointerList;
    static PerDB<T>* thePerDB;
};

#endif
