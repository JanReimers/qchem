// File: FittedCD.C  Fitted charge density: fit it to a density, then query its Coulomb matrix.
export module qchem.FittedCD;
import qchem.ChargeDensity.Types;
import qchem.ChargeDensity;            // rDM_CD (the density to fit; cross-cast to its AO face in the Imp)

export namespace qchem::ChargeDensity
{

//----------------------------------------------------------------------------------
//
//  A charge density that can be (re)fit to a density matrix and queried for its Coulomb (repulsion)
//  matrix and self energy.  The fitting itself is done by a composed FunctionFitter (hidden in
//  FittedCDImp).
//
//  NOT a ScalarFunction any more (2026-08-24).  It advertised rho_fit(r) and delivered it by forwarding to
//  the fitter, which delivered it by handing its COEFFICIENTS back to the fit basis -- the chain that got
//  deleted for exposing them.  Nothing ever called it: its sole holder, FittedVee, uses DoFit /
//  GetRepulsion / GetSelfRepulsion only, and the two forwarding lines carried a "No UT coverage" marker
//  someone had already added.  Measured before removal: zero calls across all 758 tests.
class FittedCD
{
public:
    //! "Fit me to this density."  Takes the density by its common tChargeDensity base.  EVERY finite density
    //! -- whether it carries an exact density MATRIX (IrrepCD/CompositeCD) or is a matrix-free SAD seed
    //! (NumericCD) -- presents its OWN density-fit projection via the ProjectedDensity_AO cross-cast; FittedCD
    //! is agnostic to which, and knows nothing of seeding.  See project_numericcd_refactor.
    virtual void    DoFit          (const rChargeDensity&         )      =0;
    virtual double  GetSelfRepulsion(                              ) const=0;  // 1/2 <ro(1)|1/r12|ro(2)>
    virtual rsmat_t GetRepulsion    (const odftbs_t*               ) const=0;
    // NO Clone().  It was a pure virtual whose SOLE implementation was assert(false)+nullptr -- a contract
    // clause no implementor could honour and no caller ever used (R2.7).  Its stated purpose, building a
    // polarized CD from an unpolarized one, has a REAL blocker recorded in that assert: it needs a
    // POLYMORPHIC fitter clone, or the constrained fitter gets sliced.  Restore the clause when that
    // exists -- declaring it early only made every implementor promise something none could deliver.
};

} //namespace
