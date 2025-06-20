// File: LDAVxc.H  Exact Exchange potential, only useful for plotting.
#ifndef _LDAVxc_H_
#define _LDAVxc_H_



#include "ExchangeFunctional.H"
#include "HamiltonianTerm.H"
#include <memory>

//###############################################################################
//
//  Local density exchange potential using exact charge density.
//
class LDAVxc
    : public virtual FittablePotential
    , public virtual ScalarFFClient
    , private        Dynamic_HT_Imp
{
    typedef std::shared_ptr<ExFunctional> ex_t;
public:
    LDAVxc();
    LDAVxc(ex_t& lda);
    // Required by HamiltonianTerm
    virtual void UseChargeDensity(const DM_CD* exact);
    virtual void GetEnergy       (EnergyBreakdown&,const DM_CD* cd         ) const;
    // Required by FittablePotential.
    virtual const ScalarFunction<double>* GetScalarFunction() const {return itsExchangeFunctional.get();}
    
    virtual std::ostream&           Write(std::ostream&) const;
    virtual std::istream&           Read (std::istream&)      ;
private:
    virtual SMat CalcMatrix(const ibs_t*,const Spin&,const DM_CD* cd) const;

    ex_t itsExchangeFunctional;
};

#endif //_LDAVxc_H_
