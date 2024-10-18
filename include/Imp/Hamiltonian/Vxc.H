// File: HartreeFockVxc.H  HartreeFock Exchange potential, only useful for plotting.
#ifndef _HartreeFockVxc_H_
#define _HartreeFockVxc_H_



#include "Imp/Hamiltonian/HamiltonianTerm.H"

//###############################################################################
//
//  HartreeFock exchange potential.
//
class HartreeFockVxc
    : public virtual HamiltonianTerm
    , private        HamiltonianTermImplementation
{
public:
    HartreeFockVxc();

    // Required by HamiltonianTerm
    virtual void           GetEnergy(TotalEnergy& ) const;
    virtual std::ostream&  Write    (std::ostream&) const;

private:
    virtual bool DependsOnChargeDensity() const {return true;}
    virtual SMat CalculateHamiltonianMatrix(const IrrepBasisSet*,const Spin&) const;
};

#endif //_HartreeFockVxc_H_
