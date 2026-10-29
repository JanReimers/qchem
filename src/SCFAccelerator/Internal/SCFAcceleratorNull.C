// FIle: SCFAccelerator_Null.H  A simple pass-through accelerator proxy that does no acceleration.
module;
#include <iosfwd>
#include <vector>
export module qchem.SCFAccelerator.Internal.SCFIrrepAcceleratorNull;
export import qchem.SCFAccelerator;

export namespace qchem::SCFAccelerators
{

// Plain diagonalising per-irrep accelerator (no extrapolation): transform F and SolveOrtho.
// Templated on T (rX/cX); the <double> alias preserves the existing DIIS/GDM fallback callers.
//! The scalar aggregation face of a Null block: its [F',D'] norm, whatever the block's scalar.
class NullBlock
{
public:
    virtual ~NullBlock() {}
    virtual double GetError() const=0;
};
template <class T> class tSCFIrrepAcceleratorNull : public virtual tSCFIrrepAccelerator<T>, public NullBlock
{
public:
    tSCFIrrepAcceleratorNull(const LASolver<T>* lasb,const Irrep&)
        : itsLASolver(lasb) {};
    virtual ~tSCFIrrepAcceleratorNull() {};
    //! Transform F, and measure [F',D'] (OrbitalGradient): a plain diagonalizer has no USE for the gradient, but the run
    //! does -- the convergence gate, the trace column and the linear mixer's adaptive step all read it (OpenWork §3, the
    //! [F,D] row: a Null run used to report 0, so the mixer's step froze at its first growth, alpha=1, and a symmetric
    //! 2-cycle could never be damped).
    virtual void UseFD(const hmat_t<T>& F, const hmat_t<T>& DPrime);
    virtual typename LASolver<T>::UUd_t NextOrbitals();
    virtual double GetError() const override {return itsEn;}
private:
    const LASolver<T>*   itsLASolver; //Knows the ortho transform
    hmat_t<T>            itsFPrime;
    double               itsEn=0.0;   //!< ||[F',D']|| of the last UseFD
};

//! A "no acceleration" MANAGER: each Create() hands back a plain diagonalising per-irrep Null
//! accelerator typed by the block.  Non-template per doc/RealComplexPlan.md §6 -- one manager serves
//! the real and complex paths alike.
class SCFAcceleratorNull : public virtual SCFAccelerator
{
public:
    //! The manager registers each block on its scalar aggregation face (as DIIS and GDM do) and keeps no other reference.
    virtual tSCFIrrepAccelerator<double>* Create(const LASolver<double>* las,const Irrep& qns, int)
        {auto* a=new tSCFIrrepAcceleratorNull<double>(las,qns); itsIrreps.push_back(a); return a;}
    virtual tSCFIrrepAccelerator<dcmplx>* Create(const LASolver<dcmplx>* las,const Irrep& qns, int)
        {auto* a=new tSCFIrrepAcceleratorNull<dcmplx>(las,qns); itsIrreps.push_back(a); return a;}
    virtual bool   CalculateProjections()               {return false;}
    virtual void   ShowLabels     (std::ostream&) const {}
    virtual void   ShowConvergence(std::ostream&) const {}
    //! The WORST block's [F',D'] (WorstIrrep), as every manager reports it.
    virtual double GetError       ()              const
    {
        double e=0.0;
        for (const NullBlock* k : itsIrreps) e=WorstIrrep(e, k->GetError());
        return e;
    }
    virtual const char* Tag       ()              const {return "Null";}
private:
    std::vector<const NullBlock*> itsIrreps;   //!< the scalar aggregation face -- children may differ in T
};

using SCFIrrepAcceleratorNull = tSCFIrrepAcceleratorNull<double>;

} //namespace
