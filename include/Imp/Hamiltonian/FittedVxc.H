// File: FittedVxc.H  Fitted exchange potential.
#ifndef _FittedVxc_H_
#define _FittedVxc_H_



#include "Imp/Hamiltonian/HamiltonianTerm.H"
#include "Imp/Hamiltonian/ExchangeFunctional.H"
#include "Imp/Fitting/FittedFunction.H"

//###############################################################################
//
//  Linear least squares fit the exchange potential.  The fit basis set is
//  inserted by the constructor,  and is not owned by FittedVxc, and as such
//  does not get deleted by ~FittedVee.  The LDA function is owned by Vxc.
//
class FittedVxc
    : public virtual FittedFunction
    , public virtual Dynamic_HT
    , private        Dynamic_HT_Imp
    , public         FittedFunctionImp<double>
{
    typedef Static_HT::SMat SMat;
public:
    typedef FittedFunctionImp<double>::mesh_t mesh_t;
    typedef FittedFunctionImp<double>::bs_t   bs_t;
    typedef std::shared_ptr<ExFunctional>     ex_t;

    FittedVxc();
    FittedVxc(bs_t& VxcFitBasisSet, ex_t&, mesh_t&);
    ~FittedVxc();
    // Required by HamiltonianTerm
    virtual void GetEnergy       (TotalEnergy&,const DM_CD* cd              ) const;
    // Required by FittablePotential.
    virtual void UseChargeDensity(const DM_CD* exact);

    virtual std::ostream&   Write(std::ostream&) const;
    

private:
    virtual SMat CalcMatrix(const ibs_t*,const Spin&,const DM_CD* cd) const;

    FittablePotential* itsLDAVxc; //Something to fit to.
};

#endif //_FittedVxc_H_
