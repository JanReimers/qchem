// File: Vxc.H  Hartree-Fock Exchange potential
#ifndef _Vxc_H_
#define _Vxc_H_



#include "HamiltonianTerm.H"

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
    virtual void           GetEnergy(EnergyBreakdown&,const DM_CD* cd ) const;
    virtual std::ostream&  Write    (std::ostream&) const;

private:
    virtual SMat CalcMatrix(const ibs_t*,const Spin&,const DM_CD* cd) const;
};

#endif //_Vxc_H_
