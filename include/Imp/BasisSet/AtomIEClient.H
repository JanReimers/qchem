// File: AtomIEClient.H Common IE client code for all atom basis sets and IEs.
#ifndef _AtomIEClient_H_
#define _AtomIEClient_H_

#include <IEClient.H> // <> prevents file search on this folder.
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
    
    void Init(const std::set<double>& exponents,size_t l, int m=0);
    void Init(const Vector  <double>& exponents,size_t l, int m=0);
    
    virtual size_t size() const {return es.size();}
    typedef std::tuple<int,int,int,double,double> bf_tuple;
    bf_tuple operator()(index_t i) const {return std::make_tuple(n,l,m,es(i),ns(i));}
    auto indices() const {return es.indices();}
    auto indices(size_t i) const {return es.indices(i);}
    typedef std::vector<size_t>::const_iterator const_iterator;
    
private:
    virtual double Norm(double e, size_t l) const=0;
};



#endif // _AtomIEClient_H_
