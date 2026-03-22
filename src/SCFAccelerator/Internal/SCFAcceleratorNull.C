// FIle: SCFAccelerator_Null.H  A simple pass through accerlator proxy that does no acceleration.
export module qchem.SCFAccelerator.Internal.SCFIrrepAcceleratorNull;
export import qchem.SCFAccelerator;

export class SCFIrrepAcceleratorNull : public virtual SCFIrrepAccelerator
{
public:
    SCFIrrepAcceleratorNull(const LASolver<double>* las,const LASolver_blaze<double>* lasb,const Irrep_QNs&) 
        : itsLaSolver(las)
        , itsLaSolver_blaze(lasb) {};
    virtual ~SCFIrrepAcceleratorNull() {};
    virtual void UseFD(const SMatrix<double>& F, const SMatrix<double>& DPrime);
    virtual SMatrix<double> Project();
private:
    const LASolver<double>*   itsLaSolver; //Knows the ortho transform
    const LASolver_blaze<double>*   itsLaSolver_blaze; //Knows the ortho transform
    SMatrix<double>   itsFPrime;
};
