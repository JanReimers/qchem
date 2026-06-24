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
//! accelerator.  Gives the plane-wave (dcmplx) path a working diagonaliser without DIIS/GDM (which
//! remain the real-path managers).
template <class T> class tSCFAcceleratorNull : public virtual tSCFAccelerator<T>
{
public:
    virtual tSCFIrrepAccelerator<T>* Create(const LASolver<T>* las,const Irrep& qns, int)
        {return new tSCFIrrepAcceleratorNull<T>(las,qns);}
    virtual bool   CalculateProjections()               {return false;}
    virtual void   ShowLabels     (std::ostream&) const {}
    virtual void   ShowConvergence(std::ostream&) const {}
    virtual double GetError       ()              const {return 0.0;}
};

using SCFIrrepAcceleratorNull = tSCFIrrepAcceleratorNull<double>;
using SCFAcceleratorNull      = tSCFAcceleratorNull<double>;

} //namespace
