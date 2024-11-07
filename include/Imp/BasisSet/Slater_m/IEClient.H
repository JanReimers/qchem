// File: Slater_m/IEClient.H
#ifndef _Slater_m_IEClient_H_
#define _Slater_m_IEClient_H_

#include <IEClient.H> // <> prevents file search on this Slater folder.
#include "oml/vector.h"
#include "oml/matrix.h"
#include <map>
using std::cout;
using std::endl;

namespace Slater_m
{
    
struct IEData
{
    IEData() {}; 
    IEData(size_t N) : Ns(N), Ls(N), Ms(N), es(N), ns(N) {};
    Vector<size_t> Ns;  //Angular momentum.
    Vector<size_t> Ls;  //Angular momentum.
    Vector<int   > Ms;  //Angular momentum.
    Vector<double> es; //The orbital exponents.
    Vector<double> ns; //Normalization constants
};

    
struct IrrepIEClient
    : public virtual ::IrrepIEClient
    , public IEData
{
    IrrepIEClient() {};
    IrrepIEClient(size_t N) : IEData(N) {};
    void Init(double minexp,double maxexp,size_t L, int n);
    
    virtual size_t size() const {return es.size();}
    typedef std::tuple<int,int,int,double,double> bf_tuple;
    bf_tuple operator()(index_t i) const {return std::make_tuple(Ns(i),Ls(i),Ms(i),es(i),ns(i));}
};

class Cacheable
{
public:
    virtual ~Cacheable() {};
};

class SlaterCD :public virtual Cacheable
{
public:
    SlaterCD(double eab, double ecd, size_t LMax);
    Vector<double> Coulomb_Rk(int la,int lc) const; //R_k(la,la,lc,lc);
    Vector<double> ExchangeRk(int la,int lb) const; //R_k(la,lb,la,lb);
    
//private:
    static double fk(double a, double ab, int k,int n);
//    static double D (double a, double ab, int k,int n);
//    double R(int k,int la, int lb, int lc, int ld) const;
    
    double eab, ecd;
    size_t LMax;
    Matrix<double> Iab,Icd;
};

//
//
//  Cache object based on four unsigned integer indices.
//
class Cache4
{
protected:
    ~Cache4();
    void       loop_1(size_t i1) const;
    void       loop_2(size_t i2) const;
    void       loop_3(size_t i3) const;
    virtual const Cacheable* loop_4(size_t i4) const;
    
private:
    virtual const Cacheable* Create(size_t i1,size_t i2,size_t i3,size_t i4) const=0;

    typedef std::map<size_t,const Cacheable*> cache_4; 
    typedef std::map<size_t,cache_4> cache_3; 
    typedef std::map<size_t,cache_3> cache_2; 
    typedef std::map<size_t,cache_2> cache_t; 
    
    mutable cache_t cache;
    mutable cache_2* i1_cache;
    mutable cache_3* i2_cache;
    mutable cache_4* i3_cache;
    mutable size_t i1,i2,i3,i4; //Current indexes
};

//Use the same struct for a full basis set.;
struct IEClient 
    : public virtual ::IEClient 
    , public virtual Cache4
    , public IEData
{
    IEClient() {};
    IEClient(size_t N) : IEData(N) {};
    void Append(const IrrepIEClient*);
    
    virtual size_t size() const {return es.size();}
    size_t LMax() const {return Max(Ls);}
    virtual const Cacheable* Create(size_t ia,size_t ic,size_t ib,size_t id) const;
    
    void loop_1(size_t ia) const {Cache4::loop_1(es_indices[ia-1]);}
    void loop_2(size_t ic) const {Cache4::loop_2(es_indices[ic-1]);}
    void loop_3(size_t ib) const {Cache4::loop_3(es_indices[ib-1]);}

    virtual const SlaterCD* loop_4(size_t id)  const
    {
        const Cacheable* c=Cache4::loop_4(es_indices[id-1]);
        return dynamic_cast<const SlaterCD*>(c);
    }


    
    //! Get all indices with Ls(i)=l
    const std::vector<size_t>& indices(size_t l) const;
    
    std::map<double,size_t> unique_es; //Unique exponents.
    std::vector<size_t> es_indices; //zero based index for each exponent.
    std::vector<double> unique_esv;
    std::map<size_t,std::vector<size_t> > L_indices; // <L,indices>
}; 

} //namespace Slater_m

#endif // _Slater_m_IEClient_H_
