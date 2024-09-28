#ifndef PTR_VECTOR1_H_INCLUDED
#define PTR_VECTOR1_H_INCLUDED

#include <vector>
#include <cassert>

template <class T> class optr_vector1;

template <class T> class optr_vector1<T*>
    : public std::vector<T*>
{
    typedef std::vector<T*> Base;
public:
    typedef typename Base::iterator iterator;
    
    explicit optr_vector1() : Base() {};
    explicit optr_vector1(int size) : Base(size) {};
    optr_vector1(const std::vector<T*>& v) : Base(v) {}; //copy pointers
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

private:
    optr_vector1& operator=(const optr_vector1&);
};

template <class T, class D> class dynamic_cast_iterator;

template <class T, class D> class dynamic_cast_iterator<T*,D*>
{
public:
    //typedef typename std::vector<T*>::iterator ITER;
    typedef typename optr_vector1<T*>::iterator iterator;
    typedef typename optr_vector1<T*>::const_iterator const_iterator;
    
    dynamic_cast_iterator(const optr_vector1<T*>& v) : current(v.begin()) {};
    dynamic_cast_iterator(const const_iterator& i) : current(i) {};
    T* operator++(int) {current++;return *current;}
//    D* end       () const {return dcast(endp);}
    D* operator* () const {return dcast(*current);}
    D* operator->() const {return dcast(*current);}
    friend bool operator!=(const dynamic_cast_iterator& id,const const_iterator& i)
    {
        return id.current!=i;
    }
private:
    D* dcast(T* p) const
    {
        assert(dynamic_cast<D*>(p));
        return dynamic_cast<D*>(p);
    }
    const_iterator current,endp;
};


#endif // PTR_VECTOR1_H_INCLUDED
