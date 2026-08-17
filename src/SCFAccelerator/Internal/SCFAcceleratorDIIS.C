// FIle: SCFAcceleratorDIIS.H  Direct Inversion of the Iterative Subspace (DIIS) algorithm
module;
#include <limits>
#include <deque>
#include <vector>
#include <string>
#include <iostream>
#include <complex>   // std::real (the DIIS B-matrix metric is the real part of a Frobenius inner product)
export module qchem.SCFAccelerator.Internal.SCFAcceleratorDIIS;
export import qchem.SCFAccelerator;
import qchem.Blaze;

export namespace qchem::SCFAccelerators
{

// (DIISParams moved to the public qchem.SCFAccelerator 2026-08-09 -- documented there.)


//! \brief The SCALAR AGGREGATION face one per-irrep DIIS block shows its manager (doc/RealComplexPlan.md
//! §6): every manager↔child interaction is real-scalar (the B-matrix contributions, the history depth,
//! the lockstep append/purge), so the manager needs no view of the block's scalar type T at all -- which
//! is what lets ONE non-template manager coordinate real and complex blocks in the same run.
class DIIS_Block
{
public:
    virtual ~DIIS_Block() {};
    virtual size_t GetNproj() const=0;                    //!< current history depth (lockstep across blocks)
    virtual double GetError() const=0;                    //!< this block's ||[F',D']||
    virtual double GetError(size_t i, size_t j) const=0;  //!< B-matrix contribution Re tr(E_i^H E_j)
    virtual void   Append1()=0;                           //!< bank the current (F',E) into the history
    virtual void   Purge1()=0;                            //!< drop the oldest history entry
};

// Per-irrep DIIS (Pulay): cache the (F', error=[F',D']) history; the manager solves for REAL coefficients
// c_i and we extrapolate F' = Sum c_i F'_i then diagonalize.  Templated on T (rX/cX): the error metric is
// the REAL part of the Frobenius inner product <E_i,E_j> = tr(E_i^H E_j), so the B matrix and the
// coefficients c stay real for both paths -- only F'/D'/E and the LASolver carry T.  hmat_t<double> is
// rsmat_t, conj/Re are identities, so the <double> instantiation is byte-identical to the original.
template <class T> class tSCFIrrepAcceleratorDIIS
    : public virtual tSCFIrrepAccelerator<T>
    , public virtual DIIS_Block
{
public:
    typedef std::deque< mat_t<T>> mv_t; //error matrices [F',D'] (general; anti-Hermitian)
    typedef std::deque<hmat_t<T>> sv_t; //Fock matrices (Hermitian)
    typedef std::deque<double>    dv_t; //doubles

    tSCFIrrepAcceleratorDIIS(const DIISParams&,const LASolver<T>*,const Irrep&,const rvec_t& cs);
    virtual ~tSCFIrrepAcceleratorDIIS();

    virtual void UseFD(const hmat_t<T>& F, const hmat_t<T>& DPrime);
    virtual typename LASolver<T>::UUd_t NextOrbitals();
private:
    hmat_t<T> Project(); //DIIS-extrapolated (orthonormal-basis) Fock matrix.
    // The DIIS_Block scalar face (private overrides are fine: the manager calls through DIIS_Block*).
    virtual size_t GetNproj() const {return itsEs.size();}
    virtual double GetError() const {return itsEn;}
    //! B-matrix metric: Re tr(E_i^H E_j) = Re Sum_kl conj(E_i)_kl (E_j)_kl (real & symmetric, both paths).
    virtual double GetError(size_t i, size_t j) const {return std::real(blazem::sum(blazem::conj(itsEs[i]) % itsEs[j]));}
    virtual void Append1();
    virtual void Purge1();

    DIISParams itsParams;
    Irrep  itsIrrep;
    // All of these are the most recent values
    hmat_t<T>  itsFPrime,itsDPrime;
    mat_t<T>   itsE;
    double     itsEn;  // ||[F',D']||
    // Caches for F',E,En
    sv_t itsFPrimes;
    mv_t itsEs;  //Error matrices [F',D']
    dv_t itsEns; //Errors ||E||=FNorm[F',D']

    const rvec_t& itsCs;  //Projection coefficients from the manager.

    const LASolver<T>*   itsLASolver; //Knows the ortho transform
};


//! The DIIS MANAGER -- non-template (doc/RealComplexPlan.md §6): its cross-block state (the real B
//! matrix, the real Pulay coefficients, the conditioning/stuck counters) never carries a block scalar,
//! so one manager serves real and complex blocks through the typed \c Create overloads, coordinating
//! them through the \c DIIS_Block scalar face.
class SCFAcceleratorDIIS : public virtual SCFAccelerator
{
public:
    SCFAcceleratorDIIS(const DIISParams&);
    ~SCFAcceleratorDIIS();
    virtual tSCFIrrepAccelerator<double>* Create(const LASolver<double>*,const Irrep&, int occ);
    virtual tSCFIrrepAccelerator<dcmplx>* Create(const LASolver<dcmplx>*,const Irrep&, int occ);
    virtual bool   CalculateProjections();
    virtual void   ShowLabels     (std::ostream&) const;
    virtual void   ShowConvergence(std::ostream&) const;
    virtual double GetError() const;
    virtual const char* Tag  () const {return "DIIS";}
    virtual int         Count() const {return (int)GetNProj();}   // projection depth = the accel column's number
    virtual double      MinSV() const {return itsLastSVMin;}       // conditioning of the bordered B
    // Out of steam: past FDMax but unable to extrapolate (singular/tiny B) for several steps.
    virtual bool   Exhausted() const {return itsStuckCount>=3;}

private:
    // The B-matrix bordering / conditioning / coefficient solve live in the SHARED, paper-faithful engine
    // qchem.Math.DIIS (which the density-face PulayMixer also uses -- doc/SCFStrategyPlan.md §4).  The B
    // matrix + coefficients c are REAL for both rX/cX paths.
    struct md_t{rsmat_t B;double sv;double scale;};   //scale = max_i Braw(i,i): the history's own error scale

    md_t    BuildB() const;
    rsmat_t BuildPrunedB(double svmin);
    size_t  Append1();
    size_t  Purge1();
    size_t  GetNProj() const;
    bool    HasProjection() const {return itsCs.size()>=2;}
    template <class T> tSCFIrrepAccelerator<T>* CreateT(const LASolver<T>*,const Irrep&, int occ);

    DIISParams itsParams;
    std::vector<DIIS_Block*> itsIrreps;   // the scalar aggregation face -- children may differ in T

    double itsEn=0.0;
    //! Conditioning of the last bordered B.  NaN until BuildPrunedB has run at least once -- DIIS does not
    //! extrapolate (so does not build B) until [F,D] < FDMax, and this was previously read UNINITIALISED by
    //! anything that asked early.  Surfaced 2026-08-10 the moment the trace grew an svMin column.
    double itsLastSVMin=std::numeric_limits<double>::quiet_NaN();
    rvec_t itsCs;
    std::string bailoutReason;
    int itsStuckCount=0; //consecutive past-FDMax iterations with no successful extrapolation.
    bool itsSeeded=false; //true once any irrep has had a nonzero error (past the zero-density start)
};

using SCFIrrepAcceleratorDIIS = tSCFIrrepAcceleratorDIIS<double>;

} //namespace
