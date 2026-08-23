// File: BasisSet/Lattice_3D/PlaneWaveFit_IBS.C  Plane-wave auxiliary (density + potential) fit basis.
//
// The auxiliary fit partner of PlaneWave_IBS: a plane-wave {G} basis implementing BOTH the density-fit face
// cFIT_CD_ABS and the potential-fit face cFIT_SF_ABS (like the molecular EFit_IBS), over the same tunable
// {G} grid as the orbital basis (distinct instances via CreateCDFitBasisSet / CreateVxcFitBasisSet; a denser
// grid is a future joint CD+Vxc upgrade).  It exists so the plane-wave DFT-fit paths go THROUGH the orbital
// basis's factory -- never assuming orbital==fit -- even though on the orthonormal {G} basis the projection
// IS the fit and this object is otherwise computationally inert (rho-tilde/V-tilde come from the density/term;
// the Hartree/XC contraction delegates to the orbital Orbital_DFT_IBS<dcmplx>).
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
#include <complex>
export module qchem.BasisSet.Lattice_3D.PlaneWaveFit_IBS;
export import qchem.BasisSet.Fit_IBS;                    // cFIT_CD_ABS (the density-fit face)
import qchem.BasisSet.Lattice_3D.Evaluators.PW;         // PW_Evaluator base -- INTERNAL to qcLattice_BS (not re-exported)
import qchem.BasisSet.Lattice_3D.IBS;                    // EPW_Irrep_IBS<E> (the shared evaluation tier)
import qchem.BasisSet.Internal.IrrepBasisSetImp;         // GetSymmetry/GetSymt/GetIrrep + itsSymmetry
import qchem.Symmetry;                                   // sym_t (the Bloch irrep, shared with the orbital basis)
import qchem.Symmetry.Lattice_3D.SpaceGroup;             // DirectOp {W|τ} -- the direct ops (glide τ) for the raster star-average
import qchem.Types;                                      // dcmplx
import qchem.Matrix3D;                                   // Matrix3D

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
    //! carrying the Bloch irrep \a sym.  \a directOps = the crystal DIRECT point ops \f${W|\tau}\f$ (empty {} =
    //! trivial = no symmetrization) ctor-injected here, where the grid is built, for the IBZ real-space raster
    //! star-average -- \f$\tau\f$ carries the glide/screw sublattice shift for non-symmorphic crystals.
    //! \a role labels this grid's JOB in the run report ("vxcFitGrid" / "densityFit") -- a construction-time
    //! fact only the creating factory knows (CreateVxcFitBasisSet vs CreateCDFitBasisSet), stamped here so the
    //! basis can SELF-report: every created grid announces at birth, labeled, and a grid that is never used
    //! still shows up -- deliberately, so stale grid construction is visible (user ruling 2026-08-16).
    PlaneWaveFit_IBS(const PW_Grid_Evaluator& e, const sym_t& sym, const std::string& role,
                     std::vector<Symmetry::Lattice_3D::DirectOp> directOps = {})
        : BasisSet::IrrepBasisSetImp<dcmplx>(sym)
        , PW_Grid_Evaluator(e)
        , itsDirectOps(std::move(directOps))
    { AnnounceGrid(role); }

    //! \copydoc BasisSet::FIT_SF_ABS::Overlap
    //! \f$\langle e^{iG_a\cdot r}/\sqrt\Omega\,|\,f\rangle=\sqrt\Omega\,\tilde f(G_a)\f$ over THIS basis's
    //! own \f$\{G\}\f$ -- sample on my raster, forward-FFT, gather the ball.  My functions carry the
    //! \f$1/\sqrt\Omega\f$ norm (\c PW_Evaluator::Eval), which is why the projection carries
    //! \f$\sqrt\Omega\f$ and \c OverlapDiagonal() below is exactly one.
    //!
    //! \warning TRUNCATION, and that is WHY \c OrthoNormalScalarFitter does not fit through this.  The
    //! XC assembly looks \f$\tilde v\f$ up at ORBITAL index differences \f$m_i-m_j\f$, which reach twice
    //! the orbital ball -- outside this fit ball, where the honest answer here is ZERO but the raster's
    //! \c GridCoeff wraps (aliases) and returns the resolved value.  Routing that fitter through this
    //! projection would therefore silently DELETE the bulk of \f$v_{xc}\f$, not move a last bit; the
    //! same distinction \c G_FieldEvaluator::ProjectField warns about, met here for real.  So this is the
    //! representation's honest self-projection and its fitter keeps the raster path (2026-08-23).
    vec_t<dcmplx> Overlap(const ScalarFunction<double>& f) const override
    {
        const cvec_t  Vt=PW_Grid_Evaluator::ForwardFFT(f(PW_Grid_Evaluator::Mesh().Points()));
        const double  rootOmega=std::sqrt(PW_Evaluator::Volume());
        const auto&   Gs=PW_Evaluator::Gs();
        vec_t<dcmplx> p(Gs.size());
        for (size_t a=0; a<Gs.size(); a++) p[a]=rootOmega*PW_Grid_Evaluator::GridCoeff(Vt, Gs[a]);
        return p;
    }
    //! \copydoc BasisSet::FIT_SF_ABS::Charge
    //! \f$\langle e^{iG\cdot r}/\sqrt\Omega|1\rangle=\sqrt\Omega\,\delta_{G,0}\f$ -- a cell integral kills
    //! every wave but the constant one.  (So \f$\int\sum_a c_a f_a=\sqrt\Omega\,c_0\f$: the mean.)
    vec_t<dcmplx> Charge() const override
    {
        const auto& Gs=PW_Evaluator::Gs();
        vec_t<dcmplx> I(Gs.size(), dcmplx(0.0));
        for (size_t a=0; a<Gs.size(); a++)
            if (Gs[a].x==0 && Gs[a].y==0 && Gs[a].z==0) I[a]=dcmplx(std::sqrt(PW_Evaluator::Volume()));
        return I;
    }

    //! \copydoc BasisSet::FIT_SF_ABS::OverlapDiagonal
    //! ALL ONES: plane waves over the cell are orthoNORMAL, which is the special case
    //! \c OrthoNormalFunctionFitter is built on -- it asserts this contract and then never calls this,
    //! so the vector below is materialised only if someone asks (a test, or a general fitter).
    vec_t<dcmplx> OverlapDiagonal() const override {return vec_t<dcmplx>(PW_Evaluator::size(), dcmplx(1.0));}

    //! A plane-wave {G} fit basis IS orthonormal (metric = I, the projection IS the fit) -- the single override
    //! satisfying the \c isOrtho contract for BOTH the density and potential fit faces.
    bool isOrtho() const override {return true;}

    //! IBZ real-space raster star-average (doc/GPWPlan1.md items 3 + 5): ρ_sym[x]=(1/|ops|)Σ_op ρ(W·x + τ) over
    //! the crystal DIRECT ops \f${W|\tau}\f$.  The rotation \f$W\cdot x\f$ is an exact integer voxel permutation;
    //! the glide translation \f$\tau\f$ is applied EXACTLY via the FFT shift theorem (\f$\rho(x+\tau)=\mathcal
    //! F^{-1}[e^{+iG\cdot\tau}\hat\rho]\f$) -- one FFT round-trip per distinct \f$\tau\f$, so a diamond ¼ glide is
    //! exact even on the \f$n=15\f$ grid (a grid NOT commensurate with τ; nearest-voxel/interpolation would
    //! alias/blur the density and mistune the SCF).  τ=0 for a symmorphic crystal, so this collapses to the plain
    //! permutation average.  Real-space average of a (near-)non-negative raster stays ≥0, so XC keeps its ρ_DM
    //! feed.  No-op when unfolded ({} ops).  Needs the CUBIC FFT grid (axis-mixing W).
    void Symmetrize(rvec_t& rho) const override
    {
        if (itsDirectOps.empty()) return;                                  // trivial {E}: exact no-op
        const ivec3_t N = this->FFTGrid();
        const int n = N.x;
        if (N.y!=n || N.z!=n || (size_t)n*n*n != rho.size()) return;        // guard: cubic grid == the raster
        auto at=[n](long ix,long iy,long iz)->size_t
                { auto m=[n](long i){return size_t(((i%n)+n)%n);}; return (m(ix)*n+m(iy))*n+m(iz); };

        // Exact fractional shift of the (periodic, band-limited) raster by +τ via the FFT phase.  n odd here
        // (5-smooth AutoGrid) -> no Nyquist bin, so the signed-frequency phase e^{+2πi m·τ} is the exact
        // fractional-delay filter (Hermitian -> the inverse transform is real).
        const double twoPi = 8.0*std::atan(1.0);
        auto shiftBy=[&](const rvec3_t& tau)->rvec_t
        {
            cvec_t F = this->ForwardFFT(rho);                              // physical, /Npts-normalised (raster order)
            for (int ax=0; ax<n; ++ax) for (int ay=0; ay<n; ++ay) for (int az=0; az<n; ++az)
            {
                int mx=ax<=n/2?ax:ax-n, my=ay<=n/2?ay:ay-n, mz=az<=n/2?az:az-n;   // signed frequency (-n/2..n/2)
                double ph = twoPi*(mx*tau.x + my*tau.y + mz*tau.z);
                F[(size_t(ax)*n+ay)*n+az] *= std::polar(1.0, ph);          // e^{+iG·τ}
            }
            return this->BackwardFFT(F);
        };

        // Pre-shift the raster once per DISTINCT nonzero τ (diamond: a single ¼-glide τ shared by 24 ops).
        std::vector<rvec3_t> taus; std::vector<rvec_t> shifted;
        auto srcFor=[&](const rvec3_t& tau)->const rvec_t&
        {
            if (std::fabs(tau.x)<1e-9 && std::fabs(tau.y)<1e-9 && std::fabs(tau.z)<1e-9) return rho;
            for (size_t i=0;i<taus.size();++i)
                if (std::fabs(taus[i].x-tau.x)<1e-9 && std::fabs(taus[i].y-tau.y)<1e-9 && std::fabs(taus[i].z-tau.z)<1e-9)
                    return shifted[i];
            taus.push_back(tau); shifted.push_back(shiftBy(tau)); return shifted.back();
        };

        rvec_t out(rho.size(), 0.0);
        const double w = 1.0/double(itsDirectOps.size());
        for (int ix=0; ix<n; ++ix) for (int iy=0; iy<n; ++iy) for (int iz=0; iz<n; ++iz)
        {
            double s=0.0;
            for (const auto& op : itsDirectOps)
            {
                const rvec_t& src = srcFor(op.tau);                        // ρ(·+τ), pre-shifted exactly
                rvec3_t g = op.W * rvec3_t(double(ix),double(iy),double(iz));   // W·x (exact integer permutation)
                s += src[at(std::lround(g.x),std::lround(g.y),std::lround(g.z))];
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
    std::vector<Symmetry::Lattice_3D::DirectOp> itsDirectOps;   //!< crystal direct ops {W|τ} for the IBZ raster star-average ({} = no-op)
};

} //namespace
