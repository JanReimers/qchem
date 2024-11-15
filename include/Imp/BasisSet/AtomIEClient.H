// File: AtomIEClient.H Common IE client code for all atom basis sets and IEs.
#ifndef _AtomIEClient_H_
#define _AtomIEClient_H_

#include <IEClient.H> // <> prevents file search on this folder.
#include "oml/vector.h"
#include <set>

struct AtomIEClient : public virtual ::IrrepIEClient
{
    AtomIEClient() {}; 
    AtomIEClient(size_t N) : Ns(N), Ls(N), Ms(N), es(N), ns(N) {};
    Vector<size_t> Ns;  //Principle QN.
    Vector<size_t> Ls;  //Angular momentum.
    Vector<int   > Ms;  //Angular momentum.
    Vector<double> es; //The orbital exponents.
    Vector<double> ns; //Normalization constants
    
    void Init(double minexp,double maxexp,size_t L, int m=0);
    void Init(std::set<double>& exponents,size_t L, int m=0);
    
    virtual size_t size() const {return es.size();}
    typedef std::tuple<int,int,int,double,double> bf_tuple;
    bf_tuple operator()(index_t i) const {return std::make_tuple(Ns(i),Ls(i),Ms(i),es(i),ns(i));}
    auto indices() const {return es.indices();}
private:
    virtual double Norm(double e, size_t l) const=0;
};

#endif // _AtomIEClient_H_
