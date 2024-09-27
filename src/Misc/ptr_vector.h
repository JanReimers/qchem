// File: ptr_vector.h  Experimental pointer vector, derived from and STL vector<void*>.
#ifndef _ptr_vector_h_
#define _ptr_vector_h_

// Copyright (1994-2003), Jan N. Reimers

#include "oml/imp/index_t.h"
#include "Misc/void_types.h"
#include <vector>
#include <cassert>
//
// Primary template for un-owned pointers.
//
template <class T> class ptr_vector;

//
//  Iterator.  Has different symmantics than a T* iterator, op* returns a T&
//  and op& returns a T*.
//
template<class T, class Ref, class Ptr, class Base>
struct ptr_vector_iterator : public Base
{
  typedef typename VoidType<T*>::void_type void_type;
  typedef ptr_vector_iterator<T,Ref,Ptr,Base> Self;
  typedef ptr_vector_iterator<T,T& ,T* ,Base > iterator;

  ptr_vector_iterator(                       ) : Base  (        ) {};
//  template <class B> ptr_vector_iterator(const B & x) : Base(static_cast<const Base&>(x)) {} //Not type safe.
  template <class B> ptr_vector_iterator(const B & x) : Base(x) {} //Not type safe.

//  bool operator==(const Self& x) const {return itsRep==x.itsRep;}
//  bool operator!=(const Self& x) const {return itsRep!=x.itsRep;}


  Ptr& operator& () const
  {
       const void_type& cv=Base::operator*();
       void_type& v=const_cast<void_type&>(cv);
       Ptr& p=reinterpret_cast<Ptr&>(v);
       return p;
  }
  Ref  operator* () const { return *(operator&()); }
  Ptr  operator->() const { return   operator&() ; }

  Self& operator++(   ) {Base::operator++();return *this;}
  Self  operator++(int) {Self tmp = *this;++*this;return tmp;}
  Self& operator--(   ) {Base::operator--();return *this;}
  Self  operator--(int) {Self tmp = *this;--*this;return tmp;}

//  void_type* itsRep; //Don't use this!!
 private:
  friend class ptr_vector<T*>;
//  friend class ptr_vector_iterator<T,const T& ,const T* >; //Give silly error message.
//  void** itsRep; //Don't use this!!
//  ptr_vector_iterator(void_type* x) : itsRep(x) {}; //Not type safe;
//  ptr_vector_iterator(void_type const* x)
//    : itsRep(const_cast<void_type*>(x)) {}; //Not type safe.  Kludge!?!
};

/*! \class ptr_vector<T*> ptr_vector.h Misc/ptr_vector.h
  \brief STL like \c vector container specialized for un-owned pointers.

  The class is derived from \c vector<void*> in order to reduce executable size. Only the
  most important members of \c vector are overriden. The iterator class supports:
  - \c op&, \c op*, \c op->, pre and post \c op++, pre and post \c op--
 */
template <class T> class ptr_vector<T*>
: public std::vector<typename VoidType<T*>::void_type>
{
    public:
  typedef T*    value_type;
  typedef       value_type& reference;
  typedef const value_type& const_reference;
  typedef typename VoidType<T*>::void_type void_type;
  typedef std::vector<void_type> Base;
  typedef typename Base::iterator BaseIterator;
  typedef typename Base::reverse_iterator BaseRIterator;
public:

  typedef T     element_type;

  typedef ptr_vector_iterator <T,T&      ,      T*, BaseIterator> iterator;
  typedef ptr_vector_iterator <T,const T&,T* const, typename Base::const_iterator> const_iterator;

//  typedef std::reverse_iterator<iterator> reverse_iterator;
//  typedef std::reverse_iterator<const_iterator> const_reverse_iterator;
    typedef ptr_vector_iterator <T,T&      ,      T*, BaseRIterator> reverse_iterator;
    typedef ptr_vector_iterator <T,const T&,T* const, std::reverse_iterator<typename Base::const_iterator> > const_reverse_iterator;

  explicit ptr_vector() : Base() {};
  explicit ptr_vector(int size) : Base(size) {};
  ~ptr_vector() {};  // Somebody else owns the pointers.

  //! Get an STL like read/write iterator.
  iterator       begin()       { return Base::begin();}
  //! Get an STL like read only iterator.
  const_iterator begin() const { return Base::begin();}
  //! Get an STL like read/write iterator.
  iterator       end  ()
  {
      typename Base::iterator i=Base::end();
      return i;
  }
  //! Get an STL like read only iterator.
  const_iterator end  () const
  {
       return Base::end  ();
  }
  //! Get an STL like read/write iterator.
  reverse_iterator       rbegin()       { return Base::rbegin();}
  //! Get an STL like read only iterator.
  const_reverse_iterator rbegin() const { return Base::rbegin();}
  //! Get an STL like read/write iterator.
  reverse_iterator       rend  ()       { return Base::rend  ();}
  //! Get an STL like read only iterator.
  const_reverse_iterator rend  () const { return Base::rend  ();}

  //! Is the vector empty()?
  using Base::empty;
  //! Get number of elements.
  using Base::size;
  //! resize preserving existing elements.
  using Base::resize;
  //! Clears out the vector, \b does \b not delete the objects.
  void clear() {Base::clear();}  // Somebody else owns the pointers.
  //! Clear but do not delete, one element.
  iterator erase(iterator p)
    {
      typename Base::iterator i(p);
      return Base::erase(i);
    }

  //! Insert a range elements at p1.
  void insert(iterator p1, iterator first, iterator last)
  {
      typename Base::iterator ip1   (p1   );
      typename Base::iterator ifirst(first);
      typename Base::iterator ilast (last );
      Base::insert(ip1,ifirst,ilast);
  }
   void insert(iterator p1, T* ptr)
  {
      typename Base::iterator ip1   (p1);
      Base::insert(ip1,ptr);
  }
  //! Non const pointer to first element.
  reference        front()       { return &begin(); }
  //! const pointer to last element.
  const value_type front() const { return const_cast<const value_type>(&begin()); }
  //! Non const pointer to first element.
  reference        back ()       { return &(--end()); }
  //! const pointer to last element.
  const value_type back () const {return const_cast<const value_type>(&--end());}

  //! Add one element to the end of the list.
  void push_back (const value_type x) {Base::push_back (static_cast<void_type>(x));}

  void pop_back () {Base::pop_back();}

  //! Const array indexing.
  const_reference operator[](int i) const {return *reinterpret_cast<T*const*>(&Base::operator[](i));}
  //! Non-const array indexing.
        reference operator[](int i)       {return *reinterpret_cast<T**>(&Base::operator[](i));}

    //
    //  Support (index_t i:arr) range iterators over indices
    //
    class index_iterator
    {
    public:
        index_iterator(index_t i) : current{i} {};
        index_iterator operator++(){current++;return (*this);}
        const index_t operator*() const {return current;}
              index_t operator*() {return current;}
        bool operator!=(const index_iterator& b) {return current!=b.current;}
    private:
        index_t current;
    };

    class iterator_proxy
    {
    public:
        iterator_proxy(index_t l, index_t h) : low(l), high(h) {};
        index_iterator begin() const {return low;}
        index_iterator end  () const {return high+1;}
    private:
        index_t low;
        index_t high;
    };

    iterator_proxy indices() const {return iterator_proxy(0,size()-1);}

};

//
// Primary template for un-owned pointers.
//
template <class T> class optr_vector;
//
//  Specialize for any pointer type.
//

/*! \class optr_vector\<T*\> ptr_vector.h Misc/ptr_vector.h
  \brief STL like \c vector container specialized for \b owned pointers.
  Same as \c ptr_vector but objects get deleted when \c optr_vector does.
  Copying is not allowed, because only one \c optr_vector can own the pointers.
*/
template <class T> class optr_vector<T*>
: private ptr_vector<T*>
{
  typedef ptr_vector<T*> Base;
 public:
  using Base::element_type;

  typedef typename Base::      iterator               iterator;
  typedef typename Base::const_iterator               const_iterator;
  typedef typename Base::      reverse_iterator       reverse_iterator;
  typedef typename Base::const_reverse_iterator const_reverse_iterator;

  explicit optr_vector() : Base() {};
  explicit optr_vector(int size) : Base(size) {};
  optr_vector(const ptr_vector<T*>&); //copy pointers
  optr_vector(const optr_vector&); //Needs to clone all non-null pointer.

  ~optr_vector() {clear();}
  void clear();
  void erase(iterator i);

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
  using Base::indices;
  
 private:
  optr_vector& operator=(const optr_vector&);

};

template <class T> void optr_vector<T*>::clear()
{
  for (iterator i=begin();i!=end();i++)
    if (&i) delete &i;
  Base::clear();
}

template <class T> void optr_vector<T*>::erase(iterator i)
{
  if (&i)  delete &i;
  Base::erase(i);
}

template <template <class> class vec, class T1, class T2>
  inline vec<T2> StaticCast(const vec<T1>& v)
{
  vec<T2> ret;
  for (typename vec<T1>::const_iterator i=v.begin();i!=v.end();i++) ret.push_back(static_cast<T2>(&i));
  return ret;
}

template <template <class> class vec, class T1, class T2>
  inline vec<T2> StaticCast(vec<T1>& v)
{
  vec<T2> ret;
  for (typename vec<T1>::iterator i=v.begin();i!=v.end();i++) ret.push_back(static_cast<T2>(&i));
  return ret;
}

template <template <class> class vec, class T1, class T2>
  inline vec<T2> ConstCast(const vec<T1>& v)
{
  vec<T2> ret;
  for (typename vec<T1>::const_iterator i=v.begin();i!=v.end();i++) ret.push_back(const_cast<T2>(&i));
  return ret;
}

// Type T must be clonable in order to copy construct.
template <class T> optr_vector<T*>::optr_vector(const optr_vector& other)
{
     for (const_iterator i=other.begin();i!=other.end();i++)
     {
        if (&i)
            push_back(i->Clone());
        else
            push_back(nullptr);
     }

}
template <class T> optr_vector<T*>::optr_vector(const ptr_vector<T*>& other)
{
     for (const_iterator i=other.begin();i!=other.end();i++)
            push_back(&i);

}


#endif // _ptr_vector_h_
