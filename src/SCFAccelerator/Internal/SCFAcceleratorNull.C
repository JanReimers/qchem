// FIle: SCFAccelerator_Null.H  A simple pass through accerlator proxy that does no acceleration.
export module qchem.SCFAccelerator.Internal.SCFIrrepAcceleratorNull;
export import qchem.SCFAccelerator;

export namespace qchem::SCFAccelerators
{

class SCFIrrepAcceleratorNull : public virtual SCFIrrepAccelerator
{
public:
    SCFIrrepAcceleratorNull(const LASolver<double>* lasb,const Irrep&) 
        : itsLASolver(lasb) {};
    virtual ~SCFIrrepAcceleratorNull() {};
    virtual void UseFD(const smat_t<double>& F, const smat_t<double>& DPrime);
    virtual smat_t<double> Project();
private:
    const LASolver<double>*   itsLASolver; //Knows the ortho transform
    smat_t<double>   itsFPrime;
};

} //namespace
