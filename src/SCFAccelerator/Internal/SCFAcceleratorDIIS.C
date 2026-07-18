// FIle: SCFAcceleratorDIIS.H  Direct Inversion of the Iterative Subspace (DIIS) algorithm
module;
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

struct DIISParams
{
    size_t Nproj;  //Number of terms to keep for proections.
    double EMax;   //DIIS starts when E<EMax
    double EMin;   //DIIS stops when E<EMin
    double SVTol;  //DIIS bails out when the minimum singular value of B matrix is < SVTol;
};

template <class T> class tSCFAcceleratorDIIS;

// Per-irrep DIIS (Pulay): cache the (F', error=[F',D']) history; the manager solves for REAL coefficients
// c_i and we extrapolate F' = Sum c_i F'_i then diagonalize.  Templated on T (rX/cX): the error metric is
// the REAL part of the Frobenius inner product <E_i,E_j> = tr(E_i^H E_j), so the B matrix and the
// coefficients c stay real for both paths -- only F'/D'/E and the LASolver carry T.  hmat_t<double> is
// rsmat_t, conj/Re are identities, so the <double> instantiation is byte-identical to the original.
template <class T> class tSCFIrrepAcceleratorDIIS : public virtual tSCFIrrepAccelerator<T>
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
    friend class tSCFAcceleratorDIIS<T>;
    size_t GetNproj() const {return itsEs.size();}
    double GetError() const {return itsEn;}
    //! B-matrix metric: Re tr(E_i^H E_j) = Re Sum_kl conj(E_i)_kl (E_j)_kl (real & symmetric, both paths).
    double GetError(size_t i, size_t j) const {return std::real(blazem::sum(blazem::conj(itsEs[i]) % itsEs[j]));}
    void Append1();
    void Purge1();

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


template <class T> class tSCFAcceleratorDIIS : public virtual tSCFAccelerator<T>
{
public:
    tSCFAcceleratorDIIS(const DIISParams&);
    ~tSCFAcceleratorDIIS();
    virtual tSCFIrrepAccelerator<T>* Create(const LASolver<T>*,const Irrep&, int occ);
    virtual bool   CalculateProjections();
    virtual void   ShowLabels     (std::ostream&) const;
    virtual void   ShowConvergence(std::ostream&) const;
    virtual double GetError() const;
    // Out of steam: past EMax but unable to extrapolate (singular/tiny B) for several steps.
    virtual bool   Exhausted() const {return itsStuckCount>=3;}

private:
    // The B-matrix bordering / conditioning / coefficient solve live in the SHARED, paper-faithful engine
    // qchem.Math.DIIS (which the density-face PulayMixer also uses -- doc/SCFStrategyPlan.md §4).  The B
    // matrix + coefficients c are REAL for both rX/cX paths.
    struct md_t{rsmat_t B;double sv;};

    md_t    BuildB() const;
    rsmat_t BuildPrunedB(double svmin);
    size_t  Append1();
    size_t  Purge1();
    size_t  GetNProj() const;
    bool    HasProjection() const {return itsCs.size()>=2;}

    DIISParams itsParams;
    std::vector<tSCFIrrepAcceleratorDIIS<T>*> itsIrreps;

    double itsEn,itsLastSVMin;
    rvec_t itsCs;
    std::string bailoutReason;
    int itsStuckCount=0; //consecutive past-EMax iterations with no successful extrapolation.
    bool itsSeeded=false; //true once any irrep has had a nonzero error (past the zero-density start)
};

using SCFIrrepAcceleratorDIIS = tSCFIrrepAcceleratorDIIS<double>;
using SCFAcceleratorDIIS      = tSCFAcceleratorDIIS<double>;
using cSCFAcceleratorDIIS     = tSCFAcceleratorDIIS<dcmplx>;

} //namespace
