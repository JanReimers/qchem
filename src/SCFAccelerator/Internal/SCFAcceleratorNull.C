// FIle: SCFAccelerator_Null.H  A simple pass-through accelerator proxy that does no acceleration.
module;
#include <iosfwd>
export module qchem.SCFAccelerator.Internal.SCFIrrepAcceleratorNull;
export import qchem.SCFAccelerator;

export namespace qchem::SCFAccelerators
{

// Plain diagonalising per-irrep accelerator (no extrapolation): transform F and SolveOrtho.
// Templated on T (rX/cX); the <double> alias preserves the existing DIIS/GDM fallback callers.
template <class T> class tSCFIrrepAcceleratorNull : public virtual tSCFIrrepAccelerator<T>
{
public:
    tSCFIrrepAcceleratorNull(const LASolver<T>* lasb,const Irrep&)
        : itsLASolver(lasb) {};
    virtual ~tSCFIrrepAcceleratorNull() {};
    virtual void UseFD(const hmat_t<T>& F, const hmat_t<T>& DPrime);
    virtual typename LASolver<T>::UUd_t NextOrbitals();
private:
    const LASolver<T>*   itsLASolver; //Knows the ortho transform
    hmat_t<T>            itsFPrime;
};

//! A "no acceleration" MANAGER: each Create() hands back a plain diagonalising per-irrep Null
//! accelerator typed by the block.  Non-template per doc/RealComplexPlan.md §6 -- one manager serves
//! the real and complex paths alike.
class SCFAcceleratorNull : public virtual SCFAccelerator
{
public:
    virtual tSCFIrrepAccelerator<double>* Create(const LASolver<double>* las,const Irrep& qns, int)
        {return new tSCFIrrepAcceleratorNull<double>(las,qns);}
    virtual tSCFIrrepAccelerator<dcmplx>* Create(const LASolver<dcmplx>* las,const Irrep& qns, int)
        {return new tSCFIrrepAcceleratorNull<dcmplx>(las,qns);}
    virtual bool   CalculateProjections()               {return false;}
    virtual void   ShowLabels     (std::ostream&) const {}
    virtual void   ShowConvergence(std::ostream&) const {}
    virtual double GetError       ()              const {return 0.0;}
    virtual const char* Tag       ()              const {return "Null";}
};

using SCFIrrepAcceleratorNull = tSCFIrrepAcceleratorNull<double>;

} //namespace
