// File: AtomIEClient.H Common IE client code for all atom basis sets and IEs.
#ifndef _AtomIEClient_H_
#define _AtomIEClient_H_

#include <IEClient.H> // <> prevents file search on this folder.
#include "Imp/BasisSet/BFGrouper.H"
#include "Imp/Integrals/Cache4.H"
#include "oml/vector.h"
#include <set>

struct AtomIrrepIEClient : public virtual ::IrrepIEClient
{
    AtomIrrepIEClient() {}; 
    AtomIrrepIEClient(size_t N) :  es(N), ns(N) {};
    Vector<double> es; //The orbital exponents.
    Vector<double> ns; //Normalization constants
    std::vector<size_t> es_indices; //Unique exponent index
    
    size_t n,l;
    int m;
    
    void Init(double minexp,double maxexp,size_t l, int m=0);
    void Init(const std::set<double>& exponents,size_t l, int m=0);
    void Init(const Vector<double>& exponents,size_t l, int m=0);
    
    virtual size_t size() const {return es.size();}
    typedef std::tuple<int,int,int,double,double> bf_tuple;
    bf_tuple operator()(index_t i) const {return std::make_tuple(n,l,m,es(i),ns(i));}
    auto indices() const {return es.indices();}
    typedef std::vector<size_t>::const_iterator const_iterator;
    
private:
    virtual double Norm(double e, size_t l) const=0;
};

struct AtomIEClient 
    : public virtual ::IEClient 
    , public virtual Cache4
    , public BFGrouper
{
    void Append(AtomIrrepIEClient*);
    
    virtual size_t size() const;
    void loop_1(size_t ia) const {Cache4::loop_1(ia);}
    void loop_2(size_t ic) const {Cache4::loop_2(ic);}
    void loop_3(size_t ib) const {Cache4::loop_3(ib);}
    virtual Vector<double> loop_4_direct  (size_t id, size_t la, size_t lc) const=0;
    virtual Vector<double> loop_4_exchange(size_t id, size_t la, size_t lc) const=0;
    virtual const AtomIrrepIEClient* operator[](size_t i) const {return itsIrreps[i-1];}
    
private:
    std::vector<const AtomIrrepIEClient*> itsIrreps;
}; 


#endif // _AtomIEClient_H_
