// File: Iterators.H  helps class to support range based iteration
#ifndef _Iterators_H_
#define _Iterators_H_

#include <cassert>
//
// THis iterator will dynamic cast to type D
//
template <class D, class it_t> class D_iterator
{
public:
    D_iterator(const it_t& b) : current(b) {};
    D_iterator(const it_t& b, const D* c) : current(b) //STL won't let us construct current(c).
    {
        while (this->operator*()!=c) 
            ++current; //Clunky
    };
    it_t operator++() {return ++current;} //Prefix only.
    D* operator*() const
    {
        D* ret(dynamic_cast<D*>((*current).get()));
        assert(ret);
        return ret;
    }
    friend bool operator!=(const D_iterator& a, const D_iterator& b)
    {
        return a.current!=b.current;
    }
private:
    it_t current;
};

template <class D, class cit_t> class D_iterator_proxy
{
    typedef D_iterator<D,cit_t> it_t;
public:
    D_iterator_proxy(const cit_t& b, const cit_t& e) : ib(b), ie(e) {};
    D_iterator_proxy(const cit_t& b, const cit_t& e, const D* start) : ib(b,start), ie(e) {};
    it_t begin() const {return ib;}
    it_t end  () const {return ie;}
private:
    it_t ib,ie;
};


#endif //_Iterators_H_
