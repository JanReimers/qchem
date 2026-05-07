// File: FittedCD.C  Fitted charge density.
export module qchem.FittedCD;
import qchem.ChargeDensity.Types;
import qchem.FittedFunction;

export namespace qchem::ChargeDensity
{

//----------------------------------------------------------------------------------
//
//
class FittedCD
    : public virtual ScalarFunction<double>
    , public virtual Fitting::FittedFunction
{
public:
    using FittedFunction::DoFit;
    
    virtual double  GetSelfRepulsion(               ) const=0;  // 1/2 <ro(1) | 1/r12 | ro(2)>
    virtual rsmat_t GetRepulsion    (const odftbs_t*) const=0;
    //Required for creating a polarized CD from and un-polarized CD
    virtual FittedCD*  Clone  (        ) const=0;
};

} //namespace