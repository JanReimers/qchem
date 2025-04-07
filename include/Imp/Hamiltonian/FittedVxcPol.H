// File: FittedVxcPol.H  Fitted polarized exchange potential.
#ifndef _FittedVxcPol_H_
#define _FittedVxcPol_H_



#include "Imp/Hamiltonian/ExchangeFunctional.H"
#include "Imp/Hamiltonian/Hamiltonian.H"
#include <memory>

//###############################################################################
//
//  Linear least squares fit the exchange potential.  The fit basis set is
//  inserted by the constructor,  and is not owned by PolarizedFittedVxc, and as such
//  does not get deleted by ~FittedVee.  The LDA function is owned by Vxc.
//
class FittedVxcPol
    : public virtual Dynamic_HT
    , private        Dynamic_HT_Imp
{
public:
    typedef std::shared_ptr<const Mesh>          mesh_t;
    typedef std::shared_ptr<const Fit_IBS>       bs_t;
    typedef std::shared_ptr<      ExFunctional>  ex_t;
    
    FittedVxcPol();
    FittedVxcPol(bs_t&, ex_t&, mesh_t& );
   ~FittedVxcPol();
    // Required by HamiltonianTerm
    virtual void GetEnergy       (TotalEnergy&,const DM_CD* cd         ) const;

    virtual std::ostream&   Write(std::ostream&) const;
    virtual std::istream&   Read (std::istream&)      ;

private:
    virtual SMat CalculateHamiltonianMatrix(const ibs_t*,const Spin&,const DM_CD* cd) const;

    Dynamic_HT* itsUpVxc  ; //Spin up.
    Dynamic_HT* itsDownVxc; //Spin down.

};

#endif //_FittedVxcPol_H_
