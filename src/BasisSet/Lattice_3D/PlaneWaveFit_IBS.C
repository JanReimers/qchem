// File: BasisSet/Lattice_3D/PlaneWaveFit_IBS.C  Plane-wave auxiliary (density + potential) fit basis.
//
// The auxiliary fit partner of PlaneWave_IBS: a plane-wave {G} basis implementing BOTH the density-fit face
// cFIT_CD_ABS and the potential-fit face cFIT_SF_ABS (like the molecular EFit_IBS), over the same tunable
// {G} grid as the orbital basis (distinct instances via CreateCDFitBasisSet / CreateVxcFitBasisSet; a denser
// grid is a future joint CD+Vxc upgrade).  It exists so the plane-wave DFT-fit paths go THROUGH the orbital
// basis's factory -- never assuming orbital==fit -- even though on the orthonormal {G} basis the projection
// IS the fit and this object is otherwise computationally inert (rho-tilde/V-tilde come from the density/term;
// the Hartree/XC contraction delegates to the orbital Band_FT_IBS).
//
// It carries NO grid logic of its own: it IS-A PW_Evaluator (a copy of the orbital basis's grid engine) and
// the evaluator-templated EPW_Irrep_IBS<E> mixin supplies op()/Gradient/GetNumFunctions.  Both fit faces are
// empty markers over IrrepBasisSet<dcmplx> (the Coulomb/overlap metrics live in the NonOrtho refinements,
// which an orthonormal {G} basis omits), so there are ZERO stubs here.
module;
#include <iosfwd>
#include <string>
#include <complex>   // std::complex operator* (dcmplx arithmetic in EvalField) -- not visible via imports alone
export module qchem.BasisSet.Lattice_3D.PlaneWaveFit_IBS;
export import qchem.BasisSet.Fit_IBS;                    // cFIT_CD_ABS (the density-fit face)
export import qchem.BasisSet.Lattice_3D.Evaluators.PW;   // PW_Evaluator (the shared grid engine, copied in)
import qchem.BasisSet.Lattice_3D.IBS;                    // EPW_Irrep_IBS<E> (the shared evaluation tier)
import qchem.BasisSet.Internal.IrrepBasisSetImp;         // GetSymmetry/GetSymt/GetIrrep + itsSymmetry
import qchem.BasisSet.G_FieldEvaluator;                  // G_FieldEvaluator (the DIP seam it implements via B)
import qchem.Math.GMap;                                  // ΔG_Map (the coefficient map to inverse-transform)
import qchem.Symmetry;                                   // sym_t (the Bloch irrep, shared with the orbital basis)
import qchem.Math;                                       // cos, sin
import qchem.Vector3D;                                   // rvec3_t dot (operator*), scalar*vector, +=
import qchem.Types;                                      // dcmplx

export namespace qchem::BasisSet::Lattice_3D
{

//! \brief Plane-wave auxiliary fit basis: BOTH \c cFIT_CD_ABS (density) and \c cFIT_SF_ABS (potential) over a
//! plane-wave grid, sharing the orbital basis's PW_Evaluator.  A thin shell -- all evaluation comes from the
//! evaluator via the mixin.
class PlaneWaveFit_IBS
    : public virtual BasisSet::cFIT_CD_ABS            // FIT_CD_ABS<dcmplx> : IrrepBasisSet<dcmplx> (density-fit face)
    , public virtual BasisSet::cFIT_SF_ABS            // FIT_SF_ABS<dcmplx> : IrrepBasisSet<dcmplx> (potential-fit face)
    , public virtual BasisSet::G_FieldEvaluator       // inverse-transform a ΔG_Map to real space (the DIP seam)
    , public         EPW_Irrep_IBS<PW_Evaluator>      // op()/Gradient/GetNumFunctions from the evaluator
    , public         BasisSet::IrrepBasisSetImp<dcmplx> // GetSymmetry/GetSymt/GetIrrep
    , public         PW_Evaluator                     // the grid engine (Cast() target), copied from the orbital basis
{
public:
    //! Build over a copy of the orbital basis's grid engine, carrying the same Bloch irrep \a sym.
    PlaneWaveFit_IBS(const PW_Evaluator& e, const sym_t& sym)
        : BasisSet::IrrepBasisSetImp<dcmplx>(sym)
        , PW_Evaluator(e)
    {}

    //! A plane-wave {G} fit basis IS orthonormal (metric = I, the projection IS the fit) -- the single override
    //! satisfying the \c isOrtho contract for BOTH the density and potential fit faces.
    bool isOrtho() const override {return true;}

    //! G_FieldEvaluator: inverse-transform the fitter's Hermitian coefficients c(dm) to the real field
    //! f(r)=Re Σ_dm c(dm) e^{i(B·dm)·r}, using our OWN grid engine's B (GetGCartesian) -- no reciprocal
    //! lattice leaks into the fitter.  Loops the sparse ΔG_Map directly (a point evaluation, not an FFT).
    double EvalField(const ΔG_Map& c, const rvec3_t& r) const override
    {
        dcmplx s(0.0);
        for (const auto& kv : c)
        {
            rvec3_t G = GetGCartesian(kv.first);
            double  ph = G*r;
            s += kv.second * dcmplx(cos(ph), sin(ph));
        }
        return s.real();
    }
    //! ∇f(r) = Σ_dm (B·dm)·(-Im[c(dm) e^{i(B·dm)·r}]) -- the gradient of the real inverse transform.
    rvec3_t EvalFieldGradient(const ΔG_Map& c, const rvec3_t& r) const override
    {
        rvec3_t g(0.0,0.0,0.0);
        for (const auto& kv : c)
        {
            rvec3_t G = GetGCartesian(kv.first);
            double  ph = G*r;
            dcmplx  ce = kv.second * dcmplx(cos(ph), sin(ph));
            g += (-ce.imag()) * G;
        }
        return g;
    }

    virtual std::string   Name      () const override {return "PlaneWaveFit";}
    virtual std::string   BasisSetID() const override {return Name()+PW_Evaluator::IDFragment();}
    virtual std::ostream& Write     (std::ostream& os) const override
        {return os << Name() << " fit IBS: " << PW_Evaluator::size() << " plane waves";}
};

} //namespace
