// File: SCFAcceleratorGDM.C  Geometric Direct Minimization (Grassmann-manifold SCF).
//
// First cut: preconditioned steepest descent on the Grassmann manifold, per irrep.
//   - gradient            g   = occupied-virtual block of the (pseudo-canonical) Fock matrix
//   - preconditioned step  d   = -g_ai / (eps_a - eps_i)   (diagonal orbital Hessian)
//   - geodesic step (t=1)  via the SVD of the tangent (Edelman-Arias-Smith geodesic)
// The preconditioner is the inverse diagonal Hessian, so the natural step length is t=1
// (a Newton step), which sidesteps a line search.  Conjugate-gradient + parallel
// transport can be layered on top of this later for faster convergence.
//
// TEMPLATED on T (double | dcmplx): the whole Grassmann machinery generalises to the COMPLEX
// Stiefel/Grassmann manifold by replacing every transpose with a conjugate-transpose and every
// tangent-space inner product <A,B> with its real part Re Tr(A^H B).  For T=double conj is a
// no-op, so tSCFAcceleratorGDM<double> is bit-identical to the original real code.  The complex
// instantiation is what lets a plane-wave / GPW SCF at a genuinely complex k polish past a DIIS
// limit cycle (and converge a degenerate open shell), which DIIS alone cannot.
//
// Ref: Van Voorhis & Head-Gordon, Mol. Phys. 100, 1713 (2002);
//      Edelman, Arias & Smith, SIAM J. Matrix Anal. Appl. 20, 303 (1998).
module;
#include <iosfwd>
#include <vector>
export module qchem.SCFAccelerator.Internal.SCFAcceleratorGDM;
export import qchem.SCFAccelerator;
import qchem.Blaze;

export namespace qchem::SCFAccelerators
{

struct GDMParams
{
    double EMax;          //Switch from a diagonalizing first step to GDM steps once [F,D] < EMax.
    double Trust=0.1;     //Trust radius: cap the largest geodesic rotation angle (radians) per step.
};

template <class T> class tSCFAcceleratorGDM;
template <class T> class tSCFIrrepAcceleratorGDM : public virtual tSCFIrrepAccelerator<T>
{
public:
    tSCFIrrepAcceleratorGDM(const GDMParams&,const LASolver<T>*,const Irrep&,int occ);
    virtual ~tSCFIrrepAcceleratorGDM();

    virtual void UseFD(const hmat_t<T>& F, const hmat_t<T>& DPrime);
    virtual typename LASolver<T>::UUd_t NextOrbitals();
    // Direct-minimization line-search hooks (see the base interface):
    //   ComputeStep() does gradient -> CG direction -> geodesic SVD (no orbital change), and
    //     returns false in the seed/diagonalize case (the caller should diagonalize);
    //   OrbitalsAt(t,commit) returns the orbitals at geodesic fraction t of the default step.
    // NextOrbitals() == ComputeStep() ? OrbitalsAt(default_t,true) : diagonalize.
    virtual bool ComputeStep();
    virtual typename LASolver<T>::UUd_t OrbitalsAt(double t, bool commit);
private:
    friend class tSCFAcceleratorGDM<T>;
    double GetError() const {return itsEn;}

    GDMParams               itsParams;
    const LASolver<T>*      itsLASolver;
    Irrep                   itsIrrep;
    size_t                  itsNocc;
    bool                    itsHaveC;   //Have we cached a set of orbitals yet?
    mat_t<T>                itsCp;      //Orthonormal-basis orbitals (n x n), columns = MOs.
    hmat_t<T>               itsFp;      //Orthonormal-basis Fock matrix.
    double                  itsEn;      //||[F',D']|| error for this irrep.
    bool                    itsActive;

    // Conjugate-gradient state (Polak-Ribiere + parallel transport on the manifold).
    bool     itsHavePrev=false;
    mat_t<T> itsPGprev, itsDprev; //prev preconditioned gradient & search direction (full tangents at itsYprev)
    mat_t<T> itsYprev;            //prev pseudo-canonical occupied orbitals (start of last geodesic)
    mat_t<T> itsUgeo, itsVtgeo;   //SVD factors of last geodesic tangent (H = Ugeo diag Vtgeo)
    rvec_t   itsSgeo;             //last geodesic angles (singular values * step length)
    double   itsDenomPrev=0.0;    //Re<G_old, P G_old> for the Polak-Ribiere denominator

    // Step state cached by ComputeStep() and consumed by OrbitalsAt().
    mat_t<T> itsScoccPC, itsScvirPC; //pseudo-canonical occupied/virtual orbitals
    mat_t<T> itsSU, itsSVt;          //SVD factors of the tangent H = itsSU diag(itsSs) itsSVt
    rvec_t   itsSs, itsSeo, itsSev;  //singular values; pc occupied/virtual orbital energies
    mat_t<T> itsSPGfull, itsSDfull;  //full tangents for the CG history
    double   itsSdenom=0.0;          //Re<G,PG> for this step
    double   itsStdef=1.0;           //default (diagonal quadratic-model) step length
};

template <class T> class tSCFAcceleratorGDM : public virtual tSCFAccelerator<T>
{
public:
    tSCFAcceleratorGDM(const GDMParams&);
    ~tSCFAcceleratorGDM();
    virtual tSCFIrrepAccelerator<T>* Create(const LASolver<T>*,const Irrep&, int occ);
    virtual bool   CalculateProjections() {return true;} //No global coupling: each irrep steps on its own.
    virtual void   ShowLabels     (std::ostream&) const;
    virtual void   ShowConvergence(std::ostream&) const;
    virtual double GetError() const;
    virtual bool   WantsLineSearch() const {return true;} //GDM is run by the direct-min loop.
private:
    GDMParams itsParams;
    std::vector<tSCFIrrepAcceleratorGDM<T>*> itsIrreps;
    double itsEn;
};

using SCFIrrepAcceleratorGDM  = tSCFIrrepAcceleratorGDM<double>;  using cSCFIrrepAcceleratorGDM = tSCFIrrepAcceleratorGDM<dcmplx>;
using SCFAcceleratorGDM       = tSCFAcceleratorGDM<double>;       using cSCFAcceleratorGDM      = tSCFAcceleratorGDM<dcmplx>;

} //namespace
