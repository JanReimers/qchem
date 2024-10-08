// File: FittedVxc.H  Fitted exchange potential.
#ifndef _FittedVxc_H_
#define _FittedVxc_H_



#include "Hamiltonian/ExchangeFunctional.H"
#include "Misc/rc_ptr.H"
#include "FunctionsImp/FittedFunctionImplementation.H"
#include "HamiltonianImplementation/HamiltonianTermImplementation.H"
#include "FunctionsImp/DumpScalarData.H"

//###############################################################################
//
//  Linear least squares fit the exchange potential.  The fit basis set is
//  inserted by the constructor,  and is not owned by FittedVxc, and as such
//  does not get deleted by ~FittedVee.  The LDA function is owned by Vxc.
//
class FittedVxc
    : public virtual FittedFunction
    , public virtual HamiltonianTerm
    , private        HamiltonianTermImplementation
    , private        DumpScalarData
    , public         FittedFunctionImplementation<double>
{
    typedef HamiltonianTerm::SMat SMat;
public:
    FittedVxc();
    FittedVxc(const rc_ptr<IrrepBasisSet>& VxcFitBasisSet, const rc_ptr<ExchangeFunctional>&, Mesh*);
    ~FittedVxc();
    // Required by HamiltonianTerm
    virtual void GetEnergy       (TotalEnergy&               ) const;
    // Required by FittablePotential.
    virtual void UseChargeDensity(const ChargeDensity* exact);

    virtual std::ostream&   Write(std::ostream&) const;
    virtual std::istream&   Read (std::istream&)      ;

private:
    virtual bool DependsOnChargeDensity() const {return true;}
    virtual SMat CalculateHamiltonianMatrix(const IrrepBasisSet*,const Spin&) const;

    FittablePotential* itsLDAVxc; //Something to fit to.
};

#endif //_FittedVxc_H_
