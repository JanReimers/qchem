// File: BasisSet/Lattice_3D/PlaneWave_IBS.C  Plane-wave irrep basis set for one k-point.
//
// A complex (dcmplx) Orbital_1E_IBS whose functions are the normalised plane waves
// e^{i(k+G).r}/sqrt(V) for the reciprocal lattice vectors G in the cutoff set
// { G : 1/2 |k+G|^2 < Ecut }.  The wave-vector k labels the Bloch (translational)
// symmetry of the block (BlochFactory).
//
// The pure grid geometry (the {G} set, op(r), the overlap/kinetic matrices, the reusable G-space
// potential assembly) lives in the shared PW_Evaluator; this basis IS-A PW_Evaluator, and the
// evaluator-templated EPW_Orbital1E_IBS<E> mixin forwards op()/Gradient/GetNumFunctions/MakeOverlap/
// MakeKinetic to it (the molecular EOrbital_1E_IBS<E> pattern).  The density-driven G-space assembly
// (rho-tilde -> Hartree, the FFT XC route) + the external pseudopotential stay here: they own atom/model
// data and are orbital-only (an auxiliary fit basis does not answer them).
module;
#include <functional>
#include <iosfwd>
#include <string>

export module qchem.BasisSet.Lattice_3D.PlaneWave_IBS;
export import qchem.BasisSet.Band_FT_IBS;       // the abstract G-space DFT capability (+ ΔG_Map)
export import qchem.ScalarFunction;             // ScalarFunction<double> -- the real-space oracle methods' arg
export import qchem.BasisSet.Lattice_3D.Evaluators.PW;   // PW_Evaluator (the shared grid engine / base subobject)
import qchem.BasisSet.Lattice_3D.IBS;           // EPW_Orbital1E_IBS<E> (the evaluator-templated mixins)
import qchem.BasisSet.Fit_IBS;                  // cFIT_CD_ABS (the auxiliary fit basis it creates) + qcMesh::MeshParams
import qchem.BasisSet.Internal.IrrepBasisSetImp;   // IrrepBasisSetImp<T>: GetSymmetry/GetSymt/GetIrrep
export import qchem.ReciprocalLattice;             // ctor takes a ReciprocalLattice (carries the B cell)
export import qchem.Pseudopotential.Integrals_Pseudo;    // the external-PP operator-assembly mixin (+ its models)
import qchem.Structure;
import qchem.Symmetry;                             // sym_t (the Bloch irrep handed to the ctor)
import qchem.Types;

export namespace qchem::BasisSet::Lattice_3D
{

//! \brief Plane-wave basis for a single k-point: the normalised waves
//! \f$ e^{i(k+G)\cdot r}/\sqrt V \f$ over the cutoff set \f$\{G:\tfrac12|k+G|^2<E_{cut}\}\f$.
class PlaneWave_IBS
    : public EPW_Orbital1E_IBS<PW_Evaluator>          // op()/Gradient/GetNumFunctions/MakeOverlap/MakeKinetic
    , public virtual BasisSet::Band_FT_IBS            // G-space DFT (rho-tilde -> Hartree, FFT XC)
    , public virtual Pseudopotential::Integrals_Pseudo<dcmplx> // G-space external pseudopotential assembly (V_loc + V_NL)
    , public         BasisSet::IrrepBasisSetImp<dcmplx> // supplies GetSymmetry/GetSymt/GetIrrep + itsSymmetry
    , public         PW_Evaluator                     // the shared grid engine (Cast() target for the mixins)
{
public:
    //! \brief Primary constructor: the Bloch symmetry IS the k-label (mirrors the atom IBSs, which take
    //! an abstract \c sym_t and pry out their quantum number).  The crystal momentum is read from the
    //! irrep via Symmetry::Getk; the basis owns no copy of the BZ grid.
    //! \param recip   the reciprocal lattice (its UnitCell matrix is \f$B=2\pi A^{-\top}\f$).
    //! \param irrep   the Bloch irrep (a BlochQN); \f$k\f$ = Symmetry::Lattice_3D::Getk(irrep).
    //! \param Ecut    plane-wave energy cutoff (Hartree): keep \f$G\f$ with \f$\tfrac12|k+G|^2<E_{cut}\f$.
    PlaneWave_IBS(const ReciprocalLattice& recip, const sym_t& irrep, double Ecut);

    //! \brief Convenience constructor for tests/callers that work directly in BZ-grid indices: builds
    //! the Bloch irrep \c BlochFactory(N,kIndex) and delegates to the primary constructor above.
    //! \param N       Brillouin-zone grid divisions (context for the integer k-label).
    //! \param kIndex  integer k-label; the fractional crystal momentum is \f$k = kIndex/N\f$.
    PlaneWave_IBS(const ReciprocalLattice& recip, const ivec3_t& N,
                  const ivec3_t& kIndex, double Ecut);

    // --- Band_FT_IBS capability: density-driven KS assembly in reciprocal space (orbital-only). ---
    //! \brief D-free Coulomb 3-centre tensor (delta support + diagonal kernel \f$4\pi/|G_c|^2\f$).  A density
    //! contracts \f$D\f$ against it to get \f$V_H\f$ directly.  Mirrors \c Orbital_DFT_IBS::Repulsion3C.
    virtual const G_ERI3& Repulsion3C(const BasisSet::cFIT_CD_ABS& c) const override;
    //! \brief D-free overlap 3-centre tensor (delta support, empty kernel).  Mirrors \c Orbital_DFT_IBS::Overlap3C.
    virtual const G_ERI3& Overlap3C(const BasisSet::cFIT_SF_ABS& c) const override;
    //! \brief Structure-factor assembly of a per-species radial form factor (the SAD seed density face).
    virtual ΔG_Map MakeFourierDensity(const Structure* atoms,
                          const std::function<double(int Z, double g2)>& formFactor) const override;

    // --- Real-space DFT-integration oracles (test-only): the same questions a future Band_DFT_IBS<T>
    // implementer (e.g. GPW: Gaussian orbitals, PW/FFT density) would answer, kept here as independent
    // analytic cross-checks of the shared FFT/Poisson machinery.  Production quadratures on the FIT basis's
    // grid via the G_FieldEvaluator seam; no Hamiltonian term calls these.  (The grid PRIMITIVES -- RhoOnGrid,
    // the forward transform, the <i|V|j> assembly, Integral -- now live on G_FieldEvaluator/PW_Evaluator, which
    // this basis carries; the tests reach them there.) ---
    chmat_t Overlap  (const ScalarFunction<double>& f) const;            //!< <i|f|j> (uncached; test oracle)
    chmat_t Repulsion(const ScalarFunction<double>& rho) const;          //!< <i|V_Coul[rho]|j> (test oracle)
    double  Integral (const ScalarFunction<double>& f) const;           //!< integral f d3r (test oracle)

    // --- 1E nuclear + external pseudopotential (own the atom/model data). ---
    virtual chmat_t MakeNuclear (const Structure*) const override;  //!< Bare-Coulomb structure factor.

    //! \brief Assemble any local external potential \f$\langle G|V|G'\rangle=\frac1\Omega\sum_a
    //! v(Z_a,|\Delta G|^2)e^{-i\Delta G\cdot\tau_a}\f$, \f$\Delta G\ne 0\f$ (\f$\Delta G=0\f$ dropped).
    virtual chmat_t MakeLocalPotential(const Structure* cl, const Pseudopotential::LocalPotential& loc) const override;
    //! \brief Assemble the separable (Kleinman-Bylander) nonlocal potential (rank-1 per atom, projector, m).
    virtual chmat_t MakeSeparablePotential(const Structure* cl, const Pseudopotential::SeparablePotential& v) const override;

    //! \brief Create this basis's auxiliary plane-wave density-fit basis (a distinct PlaneWaveFit_IBS over
    //! the same \f$\{G\}\f$ grid): the Band_FT_IBS factory seam a Hartree term obtains its fitter through.
    virtual BasisSet::cFIT_CD_ABS* CreateCDFitBasisSet(const Structure* cl, const qcMesh::MeshParams& mp) const override;
    //! \brief Create this basis's auxiliary plane-wave potential (Vxc) fit basis (the overlap-metric sibling):
    //! a distinct PlaneWaveFit_IBS the XC term obtains its scalar fitter through.
    virtual BasisSet::cFIT_SF_ABS* CreateVxcFitBasisSet(const Structure* cl, const qcMesh::MeshParams& mp) const override;

    virtual std::string Name      () const override {return "PlaneWave";}
    virtual std::string BasisSetID() const override; // geometry-aware cache key (Name + k, Ecut, nG)

    virtual std::ostream& Write(std::ostream&) const override;

private:
    //! Lazily-built caches of the density-free \f$\{G\}\f$ 3-centre tensors (intrinsic to \f$\{G\}\f$, so valid
    //! for the basis's lifetime regardless of the density) -- the plane-wave analogue of the \c ERI3 cache.
    //! Both share the same delta support; \c itsRepulsion3C additionally carries the diagonal Coulomb kernel.
    mutable G_ERI3 itsRepulsion3C; mutable bool itsRepBuilt=false;      //!< Coulomb (4pi/G^2 kernel)
    mutable G_ERI3 itsOverlap3C;   mutable bool itsOvlBuilt=false;      //!< overlap (empty kernel)
};

} //namespace
