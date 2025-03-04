// File: Slater_mj/IntegralEngine.H  Integral Engine for Slater function with kappa,mj QNs.
#ifndef _Slater_mj_IntegralEngine_H_
#define _Slater_mj_IntegralEngine_H_


#include "Imp/BasisSet/Slater_m/IntegralEngine.H"
#include "Imp/BasisSet/Slater_mj/IEClient.H"

class DiracIntegralTests;

namespace Slater_mj
{
    
class DiracIntegralEngine
    : public virtual AnalyticIE<double>
{
public:
    DiracIntegralEngine();
    virtual ~DiracIntegralEngine()  {};
    virtual void Append(const ::IrrepIEClient*);

    // One electron functions for DFT and HF
    virtual SMat MakeOverlap  (iec_t* a) const;
    virtual SMat MakeKinetic  (iec_t* a) const;
    virtual  Mat MakeKinetic  (iec_t* a,iec_t* b) const;
    virtual SMat MakeNuclear  (iec_t* a, const Cluster&) const;
    virtual SMat MakeRestMass (iec_t* a) const;
    virtual RVec MakeCharge   (iec_t* a) const;
    //  DFT Rho,Vac fitting functions
    virtual SMat MakeRepulsion  (iec_t* a) const; //Used for charge and Vxc fitting.
    virtual  Mat MakeRepulsion  (iec_t*,iec_t*) const; //Used for charge and Vxc fitting.
    virtual ERI3 MakeOverlap3C  (iec_t* ab,iec_t* c) const; //Used for DFT
    virtual ERI3 MakeRepulsion3C(iec_t* ab,iec_t* c) const; //Used for DFT
    // HF functions
    virtual ERI4 MakeDirect  (const ::IrrepIEClient* a, const ::IrrepIEClient* c) const;
    virtual ERI4 MakeExchange(const ::IrrepIEClient* a, const ::IrrepIEClient* b) const;
    virtual RVec Coulomb_AngularIntegrals(size_t la, size_t lc, int ma, int mc) const;
    virtual RVec ExchangeAngularIntegrals(size_t la, size_t lc, int ma, int mc) const;
    // Cach4 functions
    virtual Vector<double> loop_4_direct  (size_t id, size_t la, size_t lc) const;
    virtual Vector<double> loop_4_exchange(size_t id, size_t la, size_t lc) const;

    static SMat merge_diag(const SMat& l,const SMat& s);
    static SMat merge_off_diag(const Mat& ls);
private:
    friend class ::DiracIntegralTests;
    static const Dirac_IrrepIEClient* dcast(iec_t*);
    static RVec merge(const RVec& l,const RVec& s);
    static ERI4 merge_diag(const ERI4& LL,const ERI4& LS,const ERI4& SL,const ERI4& SS);
    static ERI4 merge_off_diag(const ERI4 &LLLL, const M4 &LSLS, const M4 &SLSL, const ERI4 &SSSS);

    AnalyticIE<double>* itsLargeIE;
    AnalyticIE<double>* itsSmallIE;
};

class Small_IntegralEngine
    : public virtual AnalyticIE<double>
    , public Slater_m::IntegralEngine
{
public: 
    virtual double Overlap  (double ea , double eb,size_t l) const;
    virtual double Nuclear(double ea, double eb,size_t l) const;
    virtual ERI4 MakeDirect  (const ::IrrepIEClient* a, const ::IrrepIEClient* c) const;
    //virtual ERI4 MakeExchange(const ::IrrepIEClient* a, const ::IrrepIEClient* b) const;
private:
    static const IrrepIEClient* dcast(iec_t*);
    ERI4 MakeDirectLS  (const IrrepIEClient* a, const IrrepIEClient* c) const;
    ERI4 MakeDirectSL  (const IrrepIEClient* a, const IrrepIEClient* c) const;
    ERI4 MakeDirectSS  (const IrrepIEClient* a, const IrrepIEClient* c) const;
public:
    M4   MakeExchangeLS(const IrrepIEClient* a, const IrrepIEClient* b) const;
    M4   MakeExchangeSL(const IrrepIEClient* a, const IrrepIEClient* b) const;
    ERI4 MakeExchangeSS(const IrrepIEClient* a, const IrrepIEClient* b) const;
};

class Small_IntegralEngine1
    : public virtual Integrals_1E<double>
    , public AtomIE_1E<double>
    // , public AtomIE_2E<double>
    , public AtomIE_RKB<double>
{
public:
    virtual  double Integral(qchem::IType,double ea , double eb,size_t l) const;
};

} //namespace

#endif //_Slater_m_IntegralEngine_H_
