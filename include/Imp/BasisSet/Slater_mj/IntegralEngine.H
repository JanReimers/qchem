// File: Slater_mj/IntegralEngine.H  Integral Engine for Slater function with kappa,mj QNs.
#ifndef _Slater_mj_IntegralEngine_H_
#define _Slater_mj_IntegralEngine_H_


#include "Imp/BasisSet/Slater_m/IntegralEngine.H"
#include "Imp/BasisSet/Slater_mj/IEClient.H"


namespace Slater_mj
{
    
class DiracIntegralEngine
    : public virtual AnalyticIE<double>
{
public:
    DiracIntegralEngine();
    // One electron functions for DFT and HF
    virtual SMat MakeOverlap  (iec_t* a) const;
    virtual SMat MakeKinetic  (iec_t* a) const;
    virtual SMat MakeNuclear  (iec_t* a, const Cluster&) const;
    virtual RVec MakeCharge   (iec_t* a) const;
    //  DFT Rho,Vac fitting functions
    virtual SMat MakeRepulsion  (iec_t* a) const; //Used for charge and Vxc fitting.
    virtual  Mat MakeRepulsion  (iec_t*,iec_t*) const; //Used for charge and Vxc fitting.
    virtual ERI3 MakeOverlap3C  (iec_t* ab,iec_t* c) const; //Used for DFT
    virtual ERI3 MakeRepulsion3C(iec_t* ab,iec_t* c) const; //Used for DFT
    // HF functions
    virtual ERI4 MakeDirect  (const IrrepIEClient* a, const IrrepIEClient* c) const;
    virtual ERI4 MakeExchange(const IrrepIEClient* a, const IrrepIEClient* b) const;
    virtual RVec Coulomb_AngularIntegrals(size_t la, size_t lc, int ma, int mc) const;
    virtual RVec ExchangeAngularIntegrals(size_t la, size_t lc, int ma, int mc) const;
    // Cach4 functions
    virtual Vector<double> loop_4_direct  (size_t id, size_t la, size_t lc) const;
    virtual Vector<double> loop_4_exchange(size_t id, size_t la, size_t lc) const;

private:
    static const Dirac_IrrepIEClient* dcast(iec_t*);
    static SMat merge_diag(const SMat& l,const SMat& s);
    static SMat merge_off_diag(const SMat& ls);
    static RVec merge(const RVec& l,const RVec& s);

    AnalyticIE<double>* itsLargeIE;
    AnalyticIE<double>* itsSmallIE;
};

class Small_IntegralEngine
    : public virtual AnalyticIE<double>
    , public Slater_m::IntegralEngine
{
public: 
    virtual double Overlap  (double ea , double eb,size_t l) const;
};

} //namespace

#endif //_Slater_m_IntegralEngine_H_
