// File: Vee.H  Electron-Electron Coulomb potential
#ifndef _Vee_H_
#define _Vee_H_

import qchem.Hamiltonian.Internal.Term;

//###############################################################################
//
//  Implementation of the Coulomb potential
//
//            /
// Vee(r_1) = | Ro(r_2)/r_12 d^3 r_2
//           /
//
// Ro is exact charge density calculated from sum(Dab*Ga*Gb) using the density
// matrix and MO basis functions.  This is the coulomb potential used in Hartree-Fock
// claculations.
//
class Vee
    : public virtual Dynamic_HT
    , private Dynamic_HT_Imp
{
public:
    Vee();
    // Required by HamiltonianTerm
    virtual void          GetEnergy(EnergyBreakdown&,const DM_CD* cd ) const;
    virtual std::ostream& Write    (std::ostream&) const;

private:
    virtual SMat CalcMatrix(const ibs_t*,const Spin&,const DM_CD* cd) const;
};

#endif //_Vee_H_
