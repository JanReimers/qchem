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
    double EMax;  //Switch from a diagonalizing first step to GDM steps once [F,D] < EMax.
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

    GDMParams               itsParams;
    const LASolver<double>* itsLASolver;
    Irrep                   itsIrrep;
    size_t                  itsNocc;
    bool                    itsHaveC;   //Have we cached a set of orbitals yet?
    mat_t<double>           itsCp;      //Orthonormal-basis orbitals (n x n), columns = MOs.
    rsmat_t                 itsFp;      //Orthonormal-basis Fock matrix.
    double                  itsEn;      //||[F',D']|| error for this irrep.
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
