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
#include <vector>
#include <cmath>
export module qchem.BasisSet.Lattice_3D.PlaneWaveFit_IBS;
export import qchem.BasisSet.Fit_IBS;                    // cFIT_CD_ABS (the density-fit face)
import qchem.BasisSet.Lattice_3D.Evaluators.PW;         // PW_Evaluator base -- INTERNAL to qcLattice_BS (not re-exported)
import qchem.BasisSet.Lattice_3D.IBS;                    // EPW_Irrep_IBS<E> (the shared evaluation tier)
import qchem.BasisSet.Internal.IrrepBasisSetImp;         // GetSymmetry/GetSymt/GetIrrep + itsSymmetry
import qchem.Symmetry;                                   // sym_t (the Bloch irrep, shared with the orbital basis)
import qchem.Types;                                      // dcmplx
import qchem.Matrix3D;                                   // Matrix3D -- the τ=0 DIRECT ops for real-space raster symmetrization

export namespace qchem::BasisSet::Lattice_3D
{

//! \brief Plane-wave auxiliary fit basis: BOTH \c cFIT_CD_ABS (density) and \c cFIT_SF_ABS (potential) over a
//! plane-wave grid, sharing the orbital basis's PW_Evaluator.  A thin shell -- all evaluation comes from the
//! evaluator via the mixin.
class PlaneWaveFit_IBS
    : public virtual BasisSet::cFIT_CD_ABS            // FIT_CD_ABS<dcmplx> : IrrepBasisSet<dcmplx> (density-fit face)
    , public virtual BasisSet::cFIT_SF_ABS            // FIT_SF_ABS<dcmplx> : IrrepBasisSet<dcmplx> (potential-fit face)
    , public         EPW_Irrep_IBS<PW_Grid_Evaluator> // op()/Gradient/GetNumFunctions from the evaluator
    , public         BasisSet::IrrepBasisSetImp<dcmplx> // GetSymmetry/GetSymt/GetIrrep
    , public         PW_Grid_Evaluator                // the DENSITY/FIT evaluator (grid + G_FieldEvaluator); Cast() target
{
public:
    //! Build over a density/fit grid evaluator (the FFT/Poisson grid the XC + SAD paths quadrature on),
    //! carrying the Bloch irrep \a sym.  \a directOps = the τ=0 DIRECT point ops \f$W\f$ (empty {} = trivial =
    //! no symmetrization) ctor-injected here, where the grid is built, for the IBZ real-space raster star-average.
    PlaneWaveFit_IBS(const PW_Grid_Evaluator& e, const sym_t& sym,
                     std::vector<Matrix3D<double>> directOps = {})
        : BasisSet::IrrepBasisSetImp<dcmplx>(sym)
        , PW_Grid_Evaluator(e)
        , itsDirectOps(std::move(directOps))
    {}

    //! A plane-wave {G} fit basis IS orthonormal (metric = I, the projection IS the fit) -- the single override
    //! satisfying the \c isOrtho contract for BOTH the density and potential fit faces.
    bool isOrtho() const override {return true;}

    //! IBZ real-space raster star-average (doc/GPWPlan1.md item 3): ρ_sym[g]=(1/|W|)Σ_W ρ[W·g mod N] over the
    //! τ=0 DIRECT ops.  Real-space average of a non-negative raster stays ≥0, so XC keeps its ρ_DM feed.  No-op
    //! when unfolded ({} ops).  Needs the CUBIC FFT grid (axis-mixing W); a non-cubic folded cell is not wired.
    void SymmetrizeRaster(rvec_t& rho) const override
    {
        if (itsDirectOps.empty()) return;                                  // trivial {E}: exact no-op
        const ivec3_t N = this->FFTGrid();
        const int n = N.x;
        if (N.y!=n || N.z!=n || (size_t)n*n*n != rho.size()) return;        // guard: cubic grid == the raster
        auto at=[n](long ix,long iy,long iz)->size_t
                { auto m=[n](long i){return size_t(((i%n)+n)%n);}; return (m(ix)*n+m(iy))*n+m(iz); };
        rvec_t out(rho.size(), 0.0);
        const double w = 1.0/double(itsDirectOps.size());
        for (int ix=0; ix<n; ++ix) for (int iy=0; iy<n; ++iy) for (int iz=0; iz<n; ++iz)
        {
            double s=0.0;
            for (const auto& W : itsDirectOps)
            {
                rvec3_t g = W * rvec3_t(double(ix),double(iy),double(iz));
                s += rho[at(std::lround(g.x),std::lround(g.y),std::lround(g.z))];
            }
            out[(size_t(ix)*n+iy)*n+iz] = w*s;
        }
        rho = std::move(out);
    }

    // G_FieldEvaluator (EvalField/EvalFieldGradient + the FFT quadrature grid engine) is inherited from
    // PW_Evaluator -- the fit basis quadratures v_xc on its OWN grid, no per-class grid logic here.

    virtual std::string   Name      () const override {return "PlaneWaveFit";}
    virtual std::string   BasisSetID() const override {return Name()+PW_Evaluator::IDFragment();}
    virtual std::ostream& Write     (std::ostream& os) const override
        {return os << Name() << " fit IBS: " << PW_Evaluator::size() << " plane waves";}

private:
    std::vector<Matrix3D<double>> itsDirectOps;   //!< τ=0 direct ops W for the IBZ raster star-average ({} = no-op)
};

} //namespace
