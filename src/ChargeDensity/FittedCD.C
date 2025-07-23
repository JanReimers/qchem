// File: FittedCD.C  Fitted charge density.
export module qchem.FittedCD;

export import qchem.DFT_IBS;
import qchem.FittedFunction;
//----------------------------------------------------------------------------------
//
//
export class FittedCD
    : public virtual ScalarFunction<double>
    , public virtual FittedFunction
{
public:
    using FittedFunction::DoFit;
    typedef SMatrix<double> SMat;
    
    virtual double GetSelfRepulsion    (                       ) const=0;  // 1/2 <ro(1) | 1/r12 | ro(2)>
    virtual SMatrix<double>   GetRepulsion(const TOrbital_DFT_IBS<double>*) const=0;
    //Required for creating a polarized CD from and un-polarized CD
    virtual FittedCD*  Clone  (        ) const=0;
};
