// File: SphericalGaussian_m/IEClient.H
#ifndef _SphericalGaussian_m_IEClient_H_
#define _SphericalGaussian_m_IEClient_H_

#include <IEClient.H> // <> prevents file search on this SphericalGaussian_m folder.
#include "Imp/Integrals/SphericalGaussianCD.H"
#include "Imp/BasisSet/BFGrouper.H"
#include "oml/vector.h"
#include <map>
#include <set>


namespace SphericalGaussian_m
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
    void Init(std::set<double>& exponents,size_t L, int n);
    
    virtual size_t size() const {return es.size();}
    typedef std::tuple<int,int,int,double,double> bf_tuple;
    bf_tuple operator()(index_t i) const {return std::make_tuple(Ns(i),Ls(i),Ms(i),es(i),ns(i));}
};


//Use the same struct for a full basis set.;
struct IEClient 
    : public virtual ::IEClient 
    , public virtual Cache4
    , public IEData
    , public BFGrouper
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

    virtual const SphericalGaussianCD* loop_4(size_t id)  const
    {
        const Cacheable* c=Cache4::loop_4(es_indices[id-1]);
        return dynamic_cast<const SphericalGaussianCD*>(c);
    }


    
}; 

} //namespace Slater_m

#endif // _Slater_m_IEClient_H_
