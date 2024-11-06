// File: Slater_m/IEClient.H
#ifndef _Slater_m_IEClient_H_
#define _Slater_m_IEClient_H_

#include <IEClient.H> // <> prevents file search on this Slater folder.
#include "oml/vector.h"
#include "oml/matrix.h"
#include <map>


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

class SlaterCD
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

//Use the same struct for a full basis set.;
struct IEClient 
    : public virtual ::IEClient 
    , public IEData
{
    IEClient() {};
    IEClient(size_t N) : IEData(N) {};
    void Append(const IrrepIEClient*);
    
    virtual size_t size() const {return es.size();}
    size_t LMax() const {return Max(Ls);}
    size_t eindex(int ia, int ib, int ic, int id) const
    {
       size_t N=es_indices.size();
       return ((es_indices[ia-1]*N+es_indices[ib-1])*N+es_indices[ic-1])*N+es_indices[id-1];
    }
    //! Get all indices with Ls(i)=l
    const std::vector<size_t>& indices(size_t l) const;
    
    std::map<double,size_t> unique_es; //Unique exponents.
    std::vector<size_t> es_indices; //One based index for each exponent.
    std::map<size_t,std::vector<size_t> > L_indices; // <L,indices>
}; 

} //namespace Slater_m

#endif // _Slater_m_IEClient_H_
