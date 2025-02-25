// File: VxcPol.H  Polarized HF exchange potential
#ifndef _VxcPol_H_
#define _VxcPol_H_



#include "Imp/Hamiltonian/ExchangeFunctional.H"
#include "Imp/Hamiltonian/HamiltonianTerm.H"

//###############################################################################
//
//  Linear least squares fit the exchange potential.  The fit basis set is
//  inserted by the constructor,  and is not owned by PolarizedHartreeFockVxc, and as such
//  does not get deleted by ~FittedVee.  The LDA function is owned by Vxc.
//
class VxcPol
    : public virtual HamiltonianTerm
    , private        HamiltonianTermImp
{
public:
    VxcPol();
    ~VxcPol();
    // Required by HamiltonianTerm
    virtual void           GetEnergy(TotalEnergy& ) const;
    virtual std::ostream&  Write    (std::ostream&) const;
    virtual bool           IsPolarized(           ) const;

private:
    virtual bool DependsOnChargeDensity() const {return true;}
    virtual SMat CalculateHamiltonianMatrix(const Orbital_IBS<double>*,const Spin&) const;

};

#endif //_VxcPol_H_
