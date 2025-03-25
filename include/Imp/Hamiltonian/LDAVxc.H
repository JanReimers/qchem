// File: LDAVxc.H  Exact Exchange potential, only useful for plotting.
#ifndef _LDAVxc_H_
#define _LDAVxc_H_



#include "Imp/Hamiltonian/ExchangeFunctional.H"
#include "Imp/Hamiltonian/HamiltonianTerm.H"
#include <memory>

//###############################################################################
//
//  Local density exchange potential using exact charge density.
//
class LDAVxc
    : public virtual FittablePotential
    , public virtual ScalarFFClient
    , private        HamiltonianTermImp
{
    typedef std::shared_ptr<ExFunctional> ex_t;
public:
    LDAVxc();
    LDAVxc(ex_t& lda);
    // Required by HamiltonianTerm
    virtual void UseChargeDensity(const Exact_CD* exact);
    virtual void GetEnergy       (TotalEnergy&               ) const;
    // Required by FittablePotential.
    virtual const ScalarFunction<double>* GetScalarFunction() const {return itsExchangeFunctional.get();}
    
    virtual std::ostream&           Write(std::ostream&) const;
    virtual std::istream&           Read (std::istream&)      ;
private:
    virtual bool DependsOnChargeDensity() const {return true;}
    virtual SMat CalculateHamiltonianMatrix(const TOrbital_IBS<double>*,const Spin&) const;

    ex_t itsExchangeFunctional;
};

#endif //_LDAVxc_H_
