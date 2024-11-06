// File: Slater_m/IntegralEngine.H  Integral Engine for Slater function with l,m QNs.
#ifndef _Slater_m_IntegralEngine_H_
#define _Slater_m_IntegralEngine_H_


#include <AnalyticIE.H>
#include "oml/vector.h"
#include <map>


namespace Slater_m
{
    
class IrrepIEClient;
class IEClient;
class SlaterCD;
    
class IntegralEngine
    : public virtual AnalyticIE<double>
{

public:
    typedef AnalyticIE<double> IE;
    typedef ScalarFunction <double> SF;
    typedef IE::Mat Mat;
    typedef IE::SMat SMat;
    typedef IE::RVec RVec;

            
    IntegralEngine();
    ~IntegralEngine();
    
    virtual size_t size() const {return -1;}

    virtual RVec MakeNormalization(iec_t* a) const;
    virtual RVec MakeCharge       (iec_t* a) const;
    
    virtual  SMat MakeOverlap  (iec_t* a         ) const;
    virtual  SMat MakeKinetic  (iec_t* a         ) const;
    virtual  SMat MakeNuclear  (iec_t* a,const Cluster&   ) const;
    virtual  SMat MakeRepulsion(iec_t* a         ) const; //Used for charge and Vxc fitting.
    virtual   Mat MakeRepulsion(iec_t*,iec_t*) const; //Used for charge and Vxc fitting.

    virtual ERI3 MakeOverlap3C  (iec_t* ab,iec_t* c) const; //Used for DFT
    virtual ERI3 MakeRepulsion3C(iec_t* ab,iec_t* c) const; //Used for DFT

    virtual void Make4C(ERI4& J, ERI4& K, const ::IEClient*) const; //Used for Hartree-Fock

    virtual void Report(std::ostream&) const;
    virtual IE* Clone() const;

private:
    static const IrrepIEClient* dcast(iec_t*);
    
    typedef std::tuple<int,int,int,double,double> bf_tuple;
    SMat MakeOverlap  (iec_t* ab,const bf_tuple& c) const;
    SMat MakeRepulsion(iec_t* ab,const bf_tuple& c) const;
    
    mutable size_t CDinserts,CDlookups;
    
    static double FourPi2;
    
    typedef std::map<size_t,const SlaterCD*> cache_1; 
    typedef std::map<size_t,cache_1> cache_2; 
    typedef std::map<size_t,cache_2> cache_3; 
    typedef std::map<size_t,cache_3> cache_4; 
    
    mutable cache_4 SL_CDcache4;
    
    static cache_3&        find(int ia, cache_4&);
    static cache_2&        find(int ib, cache_3&);
    static cache_1&        find(int ic, cache_2&);
    //static const SlaterCD* find(int id, cache_1&);
    

};

} //namespace

#endif //_Slater_m_IntegralEngine_H_
