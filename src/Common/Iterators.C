// File: Iterators.H  helps class to support range based iteration
module;
#include <cassert>
#include <cstddef>
#include <iterator>
#include <utility>
#include <type_traits>
export module Common.Iterators;

export template <class D, class it_t> class D_iterator
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
    friend inline bool operator!=(const D_iterator& a, const D_iterator& b)
    {
        return a.current!=b.current;
    }
private:
    it_t current;
};

export template <class D, class cit_t> class D_iterator_proxy
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

//----------------------------------------------------------------------------
//
//  Generic forward iterator over any container C that exposes random indexed
//  read access:
//      element_t  C::operator[](std::size_t) const     // element_t is usually a pointer
//
//  The container's storage stays completely private: this iterator only ever
//  calls operator[], so a Structure/BasisSet/... can store its elements in a
//  vector, a deque, or synthesize them on the fly.  It models
//  std::forward_iterator, so it drives range-based for AND the C++20 <ranges>
//  adaptors (views::filter, views::transform, ...).  Element access is by
//  value (a prvalue pointer), hence iterator_category caps at input.
//
export template <class C> class IndexIterator
{
public:
    using reference         = decltype(std::declval<const C&>()[std::declval<std::size_t>()]);
    using value_type        = std::remove_cvref_t<reference>;
    using difference_type   = std::ptrdiff_t;
    using iterator_concept  = std::forward_iterator_tag;
    using iterator_category = std::input_iterator_tag;

    IndexIterator() = default;
    IndexIterator(const C* c, std::size_t i) : itsC(c), itsI(i) {};

    reference      operator* () const {return (*itsC)[itsI];}
    IndexIterator& operator++()       {++itsI; return *this;}
    IndexIterator  operator++(int)    {IndexIterator t(*this); ++itsI; return t;}

    bool operator==(const IndexIterator&) const = default;
private:
    const C*    itsC=nullptr;
    std::size_t itsI=0;
};

//----------------------------------------------------------------------------
//
//  Downcasting flavour of IndexIterator: indexes C via operator[] (which
//  returns some base pointer) and yields it dynamic_cast to D*.  Same role as
//  D_iterator above, but index-based, so the container's storage type never
//  leaks into its interface -- it only has to expose operator[] and a count.
//
export template <class D, class C> class D_IndexIterator
{
public:
    using reference         = D*;
    using value_type        = D*;
    using difference_type   = std::ptrdiff_t;
    using iterator_concept  = std::forward_iterator_tag;
    using iterator_category = std::input_iterator_tag;

    D_IndexIterator() = default;
    D_IndexIterator(const C* c, std::size_t i) : itsC(c), itsI(i) {};

    D* operator*() const
    {
        D* d=dynamic_cast<D*>((*itsC)[itsI]);
        assert(d);
        return d;
    }
    D_IndexIterator& operator++()    {++itsI; return *this;}
    D_IndexIterator  operator++(int) {D_IndexIterator t(*this); ++itsI; return t;}

    bool operator==(const D_IndexIterator&) const = default;
private:
    const C*    itsC=nullptr;
    std::size_t itsI=0;
};

export template <class D, class C> class D_IndexProxy
{
public:
    D_IndexProxy(const C* c, std::size_t n) : itsC(c), itsN(n) {};
    D_IndexIterator<D,C> begin() const {return D_IndexIterator<D,C>(itsC,0);}
    D_IndexIterator<D,C> end  () const {return D_IndexIterator<D,C>(itsC,itsN);}
private:
    const C*    itsC;
    std::size_t itsN;
};

