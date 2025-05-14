// File: RestMass.H  Reast mass c^2 term for the Dirac hamiltonian.
#ifndef _RestMass_H_
#define _RestMass_H_

#include "Imp/Hamiltonian/HamiltonianTerm.H"

class RestMass
    : public virtual Static_HT
    , private        Static_HT_Imp
{
public:
    RestMass(                         );
    // Required by HamiltonianTerm
    virtual void          GetEnergy(EnergyBreakdown&,const DM_CD* cd ) const;
    virtual std::ostream& Write    (std::ostream&) const;

private:
    virtual SMat CalculateMatrix(const ibs_t*,const Spin&) const;
};

#endif //_RestMass_H_
