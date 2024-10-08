// File: LDAVxc.H  Exact Exchange potential, only useful for plotting.
#ifndef _LDAVxc_H_
#define _LDAVxc_H_



#include "Hamiltonian/ExchangeFunctional.H"
#include "HamiltonianImplementation/HamiltonianTermImplementation.H"
#include "Misc/rc_ptr.H"

//###############################################################################
//
//  Local density exchange potential using exact charge density.
//
class LDAVxc
    : public virtual FittablePotential
    , public virtual FittedFunctionClient
    , private        HamiltonianTermImplementation
{
public:
    LDAVxc();
    LDAVxc(const rc_ptr<ExchangeFunctional>& lda);
    // Required by HamiltonianTerm
    virtual void UseChargeDensity(const ChargeDensity* exact);
    virtual void GetEnergy       (TotalEnergy&               ) const;
    // Required by FittablePotential.
    virtual void   InjectOverlaps    (FittedFunction*, const IrrepBasisSet* theFitBasisSet) const;
    virtual void   InjectRepulsions  (FittedFunction*, const IrrepBasisSet* theFitBasisSet) const;
    virtual double FitGetConstraint(                                               ) const;

    virtual std::ostream&           Write(std::ostream&) const;
    virtual std::istream&           Read (std::istream&)      ;
private:
    virtual bool DependsOnChargeDensity() const {return true;}
    virtual SMat CalculateHamiltonianMatrix(const IrrepBasisSet*,const Spin&) const;

    rc_ptr<ExchangeFunctional> itsExchangeFunctional;
};

#endif //_LDAVxc_H_
