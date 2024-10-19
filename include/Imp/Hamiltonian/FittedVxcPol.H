// File: FittedVxcPol.H  Fitted polarized exchange potential.
#ifndef _FittedVxcPol_H_
#define _FittedVxcPol_H_



#include "Imp/Hamiltonian/ExchangeFunctional.H"
#include "Imp/Hamiltonian/Hamiltonian.H"
#include "Misc/rc_ptr.H"

//###############################################################################
//
//  Linear least squares fit the exchange potential.  The fit basis set is
//  inserted by the constructor,  and is not owned by PolarizedFittedVxc, and as such
//  does not get deleted by ~FittedVee.  The LDA function is owned by Vxc.
//
class FittedVxcPol
    : public virtual HamiltonianTerm
    , private        HamiltonianTermImp
{
public:
    FittedVxcPol();
    FittedVxcPol(const rc_ptr<IrrepBasisSet>&, const rc_ptr<ExFunctional>&, const rc_ptr<Mesh>& );
   ~FittedVxcPol();
    // Required by HamiltonianTerm
    virtual void UseChargeDensity(const ChargeDensity* exact);
    virtual void GetEnergy       (TotalEnergy&               ) const;
    virtual bool IsPolarized     (                           ) const;

    virtual std::ostream&   Write(std::ostream&) const;
    virtual std::istream&   Read (std::istream&)      ;

private:
    virtual bool DependsOnChargeDensity() const {return true;}
    virtual SMat CalculateHamiltonianMatrix(const IrrepBasisSet*,const Spin&) const;

    HamiltonianTerm* itsUpVxc  ; //Spin up.
    HamiltonianTerm* itsDownVxc; //Spin down.

};

#endif //_FittedVxcPol_H_
