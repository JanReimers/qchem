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
    AtomIrrepIEClient(size_t N) : Ns(N), Ls(N), Ms(N), es(N), ns(N) {};
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

struct AtomIEClient 
    : public virtual ::IEClient 
    , public virtual Cache4
    , public BFGrouper
{
    void Append(const AtomIrrepIEClient*);
    
    virtual size_t size() const {return es.size();}
    size_t LMax() const {return Max(Ls);}
    void loop_1(size_t ia) const {Cache4::loop_1(es_indices[ia-1]);}
    void loop_2(size_t ic) const {Cache4::loop_2(es_indices[ic-1]);}
    void loop_3(size_t ib) const {Cache4::loop_3(es_indices[ib-1]);}
    virtual Vector<double> loop_4_direct  (size_t id, size_t la, size_t lc) const=0;
    virtual Vector<double> loop_4_exchange(size_t id, size_t la, size_t lc) const=0;
    
    Vector<size_t> Ns;  //Angular momentum.
    Vector<size_t> Ls;  //Angular momentum.
    Vector<int   > Ms;  //Angular momentum.
    Vector<double> es; //The orbital exponents.
    Vector<double> ns; //Normalization constants
    
    virtual const AtomIrrepIEClient* operator[](size_t i) const {return itsIrreps[i-1];}
    std::vector<const AtomIrrepIEClient*> itsIrreps;
}; 


#endif // _AtomIEClient_H_
