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
// Ref: Van Voorhis & Head-Gordon, Mol. Phys. 100, 1713 (2002);
//      Edelman, Arias & Smith, SIAM J. Matrix Anal. Appl. 20, 303 (1998).
module;
#include <iosfwd>
#include <vector>
#include <blaze/Math.h>
export module qchem.SCFAccelerator.Internal.SCFAcceleratorGDM;
export import qchem.SCFAccelerator;

export namespace qchem::SCFAccelerators
{

struct GDMParams
{
    double EMax;          //Switch from a diagonalizing first step to GDM steps once [F,D] < EMax.
    double Trust=0.1;     //Trust radius: cap the largest geodesic rotation angle (radians) per step.
};

class SCFAcceleratorGDM;
class SCFIrrepAcceleratorGDM : public virtual SCFIrrepAccelerator
{
public:
    SCFIrrepAcceleratorGDM(const GDMParams&,const LASolver<double>*,const Irrep&,int occ);
    virtual ~SCFIrrepAcceleratorGDM();

    virtual void UseFD(const rsmat_t& F, const rsmat_t& DPrime);
    virtual LASolver<double>::UUd_t NextOrbitals();
private:
    friend class SCFAcceleratorGDM;
    double GetError() const {return itsEn;}

    // Split of a GDM step, for a future direct-minimization line search:
    //   ComputeStep() does the gradient -> CG direction -> geodesic SVD (no orbital change);
    //   OrbitalsAt(t,commit) returns the orbitals at geodesic fraction t of the default step.
    // NextOrbitals() == ComputeStep(); OrbitalsAt(1,true).  A driver can instead call
    // OrbitalsAt(t,false) to trial-evaluate the energy and OrbitalsAt(t*,true) to commit.
    bool ComputeStep();                                   //returns false in the diagonalize case
    LASolver<double>::UUd_t OrbitalsAt(double t, bool commit);

    GDMParams               itsParams;
    const LASolver<double>* itsLASolver;
    Irrep                   itsIrrep;
    size_t                  itsNocc;
    bool                    itsHaveC;   //Have we cached a set of orbitals yet?
     rmat_t                 itsCp;      //Orthonormal-basis orbitals (n x n), columns = MOs.
    rsmat_t                 itsFp;      //Orthonormal-basis Fock matrix.
    double                  itsEn;      //||[F',D']|| error for this irrep.
    bool                    itsActive;

    // Conjugate-gradient state (Polak-Ribiere + parallel transport on the manifold).
    bool   itsHavePrev=false;
    rmat_t itsPGprev, itsDprev; //prev preconditioned gradient & search direction (full tangents at itsYprev)
    rmat_t itsYprev;            //prev pseudo-canonical occupied orbitals (start of last geodesic)
    rmat_t itsUgeo, itsVtgeo;   //SVD factors of last geodesic tangent (H = Ugeo diag Vtgeo)
    rvec_t itsSgeo;             //last geodesic angles (singular values * step length)
    double itsDenomPrev=0.0;    //<G_old, P G_old> for the Polak-Ribiere denominator

    // Step state cached by ComputeStep() and consumed by OrbitalsAt().
    rmat_t itsScoccPC, itsScvirPC; //pseudo-canonical occupied/virtual orbitals
    rmat_t itsSU, itsSVt;          //SVD factors of the tangent H = itsSU diag(itsSs) itsSVt
    rvec_t itsSs, itsSeo, itsSev;  //singular values; pc occupied/virtual orbital energies
    rmat_t itsSPGfull, itsSDfull;  //full tangents for the CG history
    double itsSdenom=0.0;          //<G,PG> for this step
    double itsStdef=1.0;           //default (diagonal quadratic-model) step length
};

class SCFAcceleratorGDM : public virtual SCFAccelerator
{
public:
    SCFAcceleratorGDM(const GDMParams&);
    ~SCFAcceleratorGDM();
    virtual SCFIrrepAccelerator* Create(const LASolver<double>*,const Irrep&, int occ);
    virtual bool   CalculateProjections() {return true;} //No global coupling: each irrep steps on its own.
    virtual void   ShowLabels     (std::ostream&) const;
    virtual void   ShowConvergence(std::ostream&) const;
    virtual double GetError() const;
private:
    GDMParams itsParams;
    std::vector<SCFIrrepAcceleratorGDM*> itsIrreps;
    double itsEn;
};

} //namespace
