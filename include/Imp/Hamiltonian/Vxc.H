// File: Vxc.H  Hartree-Fock Exchange potential
#ifndef _Vxc_H_
#define _Vxc_H_



#include "Imp/Hamiltonian/HamiltonianTerm.H"

//###############################################################################
//
//  Hartree-Fock exchange potential.
//
class Vxc
    : public virtual HamiltonianTerm
    , private        HamiltonianTermImp
{
public:
    Vxc();

    // Required by HamiltonianTerm
    virtual void           GetEnergy(TotalEnergy& ) const;
    virtual std::ostream&  Write    (std::ostream&) const;

private:
    virtual bool DependsOnChargeDensity() const {return true;}
    virtual SMat CalculateHamiltonianMatrix(const IrrepBasisSet*,const Spin&) const;
};

#endif //_Vxc_H_
