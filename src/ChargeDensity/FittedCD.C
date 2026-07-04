// File: FittedCD.C  Fitted charge density (a ScalarFunction you can fit to a density and query).
export module qchem.FittedCD;
import qchem.ChargeDensity.Types;
import qchem.ChargeDensity;            // rDM_CD (the density to fit; cross-cast to its AO face in the Imp)
import qchem.ScalarFunction;           // ScalarFunction<double>

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
    //! "Fit me to this density."  Takes the density by its common tChargeDensity base.  A density that
    //! carries an exact AO projection (a real density MATRIX) is cross-cast to its ProjectedDensity_AO face;
    //! a pure real-space seed (rho(r) + charge, no matrix -- the SAD NumericCD) falls through to the
    //! ScalarFunction overload below, which overlap-fits it.  See project_numericcd_refactor.
    virtual void    DoFit          (const rChargeDensity&         )      =0;
    //! "Fit me to this plain real-space density rho(r), sampled on the mesh, carrying total \a charge."
    //! The seed path: a numeric (SAD) density has no density matrix, so the fit overlap-projects rho(r)
    //! onto the fit basis here -- the work the old NumericCD::GetRepulsion3C did on demand, now relocated.
    virtual void    DoFit          (const ScalarFunction<double>&, double charge) =0;
    virtual double  GetSelfRepulsion(                              ) const=0;  // 1/2 <ro(1)|1/r12|ro(2)>
    virtual rsmat_t GetRepulsion    (const odftbs_t*               ) const=0;
    //Required for creating a polarized CD from an un-polarized CD
    virtual FittedCD*  Clone  (        ) const=0;
};

} //namespace
