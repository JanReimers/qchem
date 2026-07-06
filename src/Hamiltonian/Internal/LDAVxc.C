// File: LDAVxc.C  Exact Exchange potential, only useful for plotting.
module;
#include <memory>
export module qchem.Hamiltonian.Internal.LDAVxc;
import qchem.Hamiltonian.Internal.Term;
import qchem.Hamiltonian.Internal.ExFunctional;
import qchem.Fitting.FunctionFitter;   // Fitting::ProjectedScalar_AO (the fit-callback role lives on this impl)
import qchem.ScalarFunction;
import qchem.Hamiltonian.Types;

export namespace qchem::Hamiltonian
{

//###############################################################################
//
//  Local density exchange potential using exact charge density.
//
class LDAVxc
    : public virtual rFittablePotential
    , public virtual Fitting::ProjectedScalar_AO
    , private        rDynamic_HT_Imp
{
    typedef std::shared_ptr<ExFunctional> ex_t;
public:
    LDAVxc();
    LDAVxc(ex_t& lda);
    // Required by HamiltonianTerm
    virtual void UseChargeDensity(const rChargeDensity* exact);
    virtual void GetEnergy       (EnergyBreakdown&,const rDM_CD* cd         ) const;
    // Required by rFittablePotential.
    virtual const ScalarFunction<double>* GetScalarFunction() const {return itsExchangeFunctional.get();}
    
    virtual std::ostream&           Write(std::ostream&) const;
private:
    virtual rsmat_t CalcMatrix(const robs_t*,const Spin&,const rChargeDensity* cd) const;

    ex_t itsExchangeFunctional;
};

} //namespace
