// File: Kinetic.H  Kinetic energy term for the hamiltonian.
#ifndef _Kinetic_H_
#define _Kinetic_H_

#include "Imp/Hamiltonian/HamiltonianTerm.H"

class Kinetic
    : public virtual Static_HT
    , private        Static_HT_Imp
{
public:
    Kinetic(                         );
    // Required by HamiltonianTerm
    virtual void          GetEnergy(EnergyBreakdown&,const DM_CD* cd ) const;
    virtual std::ostream& Write    (std::ostream&) const;

private:
    virtual SMat CalculateMatrix(const ibs_t*,const Spin&) const;
};

#endif //_Kinetic_H_
