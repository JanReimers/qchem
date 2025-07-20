// File: VxcPol.H  Polarized HF exchange potential
#ifndef _VxcPol_H_
#define _VxcPol_H_



#include "HamiltonianTerm.H"
import qchem.Hamiltonian.Internal.ExFunctional;

//###############################################################################
//
//  Linear least squares fit the exchange potential.  The fit basis set is
//  inserted by the constructor,  and is not owned by PolarizedHartreeFockVxc, and as such
//  does not get deleted by ~FittedVee.  The LDA function is owned by Vxc.
//
class VxcPol
    : public virtual Dynamic_HT
    , private        Dynamic_HT_Imp_NoCache
{
public:
    VxcPol();
    ~VxcPol();
    // Required by HamiltonianTerm
    virtual void           GetEnergy(EnergyBreakdown&,const DM_CD* cd ) const;
    virtual bool           IsPolarized() const {return true;}
    virtual std::ostream&  Write    (std::ostream&) const;

private:
    virtual SMat CalcMatrix(const ibs_t*,const Spin&,const DM_CD* cd) const;

};

#endif //_VxcPol_H_
