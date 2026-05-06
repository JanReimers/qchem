// File: LDAVxc.C  Exact Exchange potential, only useful for plotting.
module;
#include <memory>
export module qchem.Hamiltonian.Internal.LDAVxc;
import qchem.Hamiltonian.Internal.Term;
import qchem.Hamiltonian.Internal.ExFunctional;

export namespace qchem::Hamiltonian
{

//###############################################################################
//
//  Local density exchange potential using exact charge density.
//
class LDAVxc
    : public virtual FittablePotential
    , public virtual Fitting::ScalarFFClient
    , private        Dynamic_HT_Imp
{
    typedef std::shared_ptr<ExFunctional> ex_t;
public:
    LDAVxc();
    LDAVxc(ex_t& lda);
    // Required by HamiltonianTerm
    virtual void UseChargeDensity(const DM_CD* exact);
    virtual void GetEnergy       (EnergyBreakdown&,const DM_CD* cd         ) const;
    // Required by FittablePotential.
    virtual const ScalarFunction<double>* GetScalarFunction() const {return itsExchangeFunctional.get();}
    
    virtual std::ostream&           Write(std::ostream&) const;
private:
    virtual rsmat_t CalcMatrix(const ibs_t*,const Spin&,const DM_CD* cd) const;

    ex_t itsExchangeFunctional;
};

} //namespace
