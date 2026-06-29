// File: LDAVxc.C  Exact Exchange potential, only useful for plotting.
module;
#include <memory>
export module qchem.Hamiltonian.Internal.LDAVxc;
import qchem.Hamiltonian.Internal.Term;
import qchem.Hamiltonian.Internal.ExFunctional;
import qchem.Fitting.FunctionFitter;   // Fitting::ScalarFFClient (the fit-callback role lives on this impl)
import qchem.ScalarFunction;
import qchem.Hamiltonian.Types;

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
    virtual void UseChargeDensity(const rChargeDensity* exact);
    virtual void GetEnergy       (EnergyBreakdown&,const DM_CD* cd         ) const;
    // Required by FittablePotential.
    virtual const ScalarFunction<double>* GetScalarFunction() const {return itsExchangeFunctional.get();}
    
    virtual std::ostream&           Write(std::ostream&) const;
private:
    virtual rsmat_t CalcMatrix(const obs_t*,const Spin&,const rChargeDensity* cd) const;

    ex_t itsExchangeFunctional;
};

} //namespace
