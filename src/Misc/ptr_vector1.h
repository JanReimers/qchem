#ifndef PTR_VECTOR1_H_INCLUDED
#define PTR_VECTOR1_H_INCLUDED

//#include "oml/imp/index_t.h"
#include <vector>
//#include <cassert>
template <class T> class optr_vector1;

template <class T> class optr_vector1<T*>
    : public std::vector<T*>
{
    typedef std::vector<T*> Base;
public:
    typedef typename Base::iterator iterator;
    
    explicit optr_vector1() : Base() {};
    explicit optr_vector1(int size) : Base(size) {};
    optr_vector1(const std::vector<T*>&); //copy pointers
    optr_vector1(const optr_vector1& v1) //Needs to clone all non-null pointer.
    {
        for (auto i:v1) push_back(i->Clone());
    }

    ~optr_vector1()
    {
        clear();
    }
    
    using Base::begin;
    using Base::end;
    using Base::rbegin;
    using Base::rend;
    using Base::empty;
    using Base::size;
    using Base::front;
    using Base::back;
    using Base::push_back;
    using Base::operator[];
    using Base::insert;

    void clear()
    {
        for (auto i:*this) if (i) delete i;
        Base::clear();
    }

    void erase(iterator i)
    {
        if (*i)  delete *i;
        Base::erase(i);
    }
private:
    optr_vector1& operator=(const optr_vector1&);
};


#endif // PTR_VECTOR1_H_INCLUDED
