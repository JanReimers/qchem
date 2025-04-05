// File: Vxc.H  Hartree-Fock Exchange potential
#ifndef _Vxc_H_
#define _Vxc_H_



#include "Imp/Hamiltonian/HamiltonianTerm.H"

//###############################################################################
//
//  Hartree-Fock exchange potential.
//
class Vxc
    : public virtual Dynamic_HT
    , private        Dynamic_HT_Imp
{
public:
    Vxc();

    // Required by HamiltonianTerm
    virtual void           GetEnergy(TotalEnergy&,const Exact_CD* cd ) const;
    virtual std::ostream&  Write    (std::ostream&) const;

private:
    virtual bool DependsOnChargeDensity() const {return true;}
    virtual SMat CalculateHamiltonianMatrix(const TOrbital_IBS<double>*,const Spin&) const;
};

#endif //_Vxc_H_
