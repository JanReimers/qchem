// File: PerDB.H  Manage a global list of ID object *'s.
#ifndef _PerDB_H_
#define _PerDB_H_



#include "Misc/UniqueID/UniqueID.H"
#include "Misc/ptr_vector.h"

template <class T> class PerRef;
template <class T> class PerList;
template <class T> class PerDBBrowser;

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

    friend class PerDBBrowser<T>;
    friend class PerRef <T>;
    friend class PerList<T>;

    optr_vector<T*>  itsPointerList;
    static PerDB<T>* thePerDB;
};

template <class T> class PerDBBrowser
{
    PerDBBrowser(const PerDB<T>& bs)
        : begin(bs.itsPointerList.begin())
        , current(begin)
        , end(bs.itsPointerList.end())
    {}

    const T* operator++(int)
    {
        current++;
        return &current;
    }
    const T& operator* () const
    {
        return *current;
    }
    const T* operator->() const
    {
        return &current;
    }
    T*const & operator& ()
    {
        return &current;
    }
//  operator const BasisFunction*  () const;
    operator bool() const
    {
        return current!=end;
    }
private:
    typedef typename optr_vector<T*>::const_iterator ITER;
    ITER begin;
    ITER current;
    ITER end;

};

#endif
