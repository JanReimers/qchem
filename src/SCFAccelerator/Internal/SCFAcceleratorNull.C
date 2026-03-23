// FIle: SCFAccelerator_Null.H  A simple pass through accerlator proxy that does no acceleration.
export module qchem.SCFAccelerator.Internal.SCFIrrepAcceleratorNull;
export import qchem.SCFAccelerator;

export class SCFIrrepAcceleratorNull : public virtual SCFIrrepAccelerator
{
public:
    SCFIrrepAcceleratorNull(const LASolver_blaze<double>* lasb,const Irrep_QNs&) 
        : itsLaSolver_blaze(lasb) {};
    virtual ~SCFIrrepAcceleratorNull() {};
    virtual void UseFD(const smat_t<double>& F, const smat_t<double>& DPrime);
    virtual smat_t<double> Project();
private:
    const LASolver_blaze<double>*   itsLaSolver_blaze; //Knows the ortho transform
    smat_t<double>   itsFPrime;
};
