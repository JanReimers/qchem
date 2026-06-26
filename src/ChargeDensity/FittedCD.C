// File: FittedCD.C  Fitted charge density (a ScalarFunction you can fit to a density and query).
export module qchem.FittedCD;
import qchem.ChargeDensity.Types;
import qchem.Fitting.FunctionFitter;   // Fitting::ProjectedDensity_AO (the fit request) + ScalarFunction

export namespace qchem::ChargeDensity
{

//----------------------------------------------------------------------------------
//
//  A real-space charge density rho(r) that can be (re)fit to a density matrix and queried for its
//  Coulomb (repulsion) matrix and self energy.  The fitting itself is done by a composed FunctionFitter
//  (hidden in FittedCDImp).
//
class FittedCD
    : public virtual ScalarFunction<double>
{
public:
    virtual void    DoFit          (const Fitting::ProjectedDensity_AO&)      =0;  //!< "fit me to this density"
    virtual double  GetSelfRepulsion(                              ) const=0;  // 1/2 <ro(1)|1/r12|ro(2)>
    virtual rsmat_t GetRepulsion    (const odftbs_t*               ) const=0;
    //Required for creating a polarized CD from an un-polarized CD
    virtual FittedCD*  Clone  (        ) const=0;
};

} //namespace
