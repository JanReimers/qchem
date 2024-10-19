// File: rc_ptr.H
#ifndef Ref_Count_Pointer_H_
#define Ref_Count_Pointer_H_

// A reference-counted pointer.  The object it points to is destroyed
//  when the last reference-counted pointer goes out of scope.  This
//  is not part of the standard, but probably should be.

template<class X> class rc_ptr
{
public:
    rc_ptr(X* p = 0);
    rc_ptr(const rc_ptr<X>&);
    rc_ptr<X>& operator= (const rc_ptr<X>&);
    ~rc_ptr();

    X& operator*() const;
    X* operator->() const;

    X* get() const;               // This is a dangerous operation.
    X* reset(X* p = 0);

private:
    X* ptr;
    int* count;
};

// Definitions of inline member functions for rc_ptr<X>.  This file is
//  included by rc_ptr.h, and should never be used on its own.

template<class X>
inline rc_ptr<X>::rc_ptr(X* p)
    : ptr(p),
    count(new int(1))
{}

template<class X>
inline rc_ptr<X>::rc_ptr(const rc_ptr<X>& p)
    : ptr(p.ptr),
    count(p.count)
{
    ++(*count);
}

template<class X>
inline rc_ptr<X>& rc_ptr<X>::operator= (const rc_ptr<X>& p)
{
    if (&p != this && p.ptr != this->ptr)
    {
        if (--(*count) == 0)
        {
            delete ptr;
            delete count;
        }

        ptr = p.ptr;
        count = p.count;
        ++(*count);
    }

    return *this;
}

template<class X>
inline rc_ptr<X>::~rc_ptr()
{
    if (--(*count) == 0)
    {
        delete ptr;
        delete count;
    }
}

template<class X>
inline X& rc_ptr<X>::operator*() const
{
    return *ptr;
}

template<class X>
inline X* rc_ptr<X>::operator->() const
{
    return ptr;
}

template<class X>
inline X* rc_ptr<X>::get() const
{
    return ptr;
}

template<class X>
inline X* rc_ptr<X>::reset(X* p)
{
    if (p != ptr)
    {
        int* new_count = new int(1);
        if (--(*count) == 0)
        {
            delete ptr;
            delete count;
        }
        ptr = p;
        count = new_count;
    }

    return ptr;
}

#endif
