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
export import qchem.BasisSet.Band_DFT_IBS;     // the abstract real-space DFT-integration capability
export import qchem.BasisSet.Band_FT_IBS;       // the abstract G-space DFT capability (+ FourierMap)
export import qchem.BasisSet.Lattice_3D.Evaluators.PW;   // PW_Evaluator (the shared grid engine / base subobject)
import qchem.BasisSet.Lattice_3D.IBS;           // EPW_Orbital1E_IBS<E> (the evaluator-templated mixins)
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
    , public virtual BasisSet::Band_DFT_IBS<dcmplx>   // real-space DFT-integration (Hartree/XC)
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
    //! \brief Density Fourier coefficients \f$\tilde\rho(\Delta m)=\frac1\Omega\sum_{G_i-G_j=\Delta m}D_{ij}\f$
    //! for a density matrix \a D in THIS plane-wave block.
    virtual FourierMap MakeFourierDensity(const chmat_t& D) const override;
    //! \brief Structure-factor assembly of a per-species radial form factor (the SAD seed density face).
    virtual FourierMap MakeFourierDensity(const Structure* atoms,
                          const std::function<double(int Z, double g2)>& formFactor) const override;
    //! \brief Hartree matrix + energy directly from the density's G-space coefficients \a rho
    //! (\f$V_H=4\pi\tilde\rho/|B\Delta m|^2\f$, \f$E_H=\tfrac\Omega2\sum 4\pi|\tilde\rho|^2/G^2\f$).
    virtual chmat_t Repulsion(const FourierMap& rho, double& Eh) const override;

    // XC route (basis owns the FFTs; see Band_FT_IBS).
    virtual rvec_t     RhoOnGrid   (const FourierMap& rho) const override;
    virtual FourierMap ForwardGrid (const rvec_t& gridValues) const override; //!< forward FFT -> Vtilde(dm)
    virtual chmat_t    Overlap     (const FourierMap& Vtilde) const override;  //!< <i|V|j>=Vtilde(dm)
    virtual chmat_t    Overlap     (const rvec_t& Vgrid)  const override;      //!< = Overlap(ForwardGrid(Vgrid))
    virtual double     Integral    (const rvec_t& fgrid)   const override;

    // --- Band_DFT_IBS capability: real-space DFT integration on the basis's own uniform grid. ---
    virtual chmat_t Overlap  (const ScalarFunction<double>& f) const override;            //!< <i|f|j> (uncached)
    virtual chmat_t Repulsion(const ScalarFunction<double>& rho, double& Eh) const override; //!< <i|V_Coul[rho]|j> (uncached)
    virtual double  Integral         (const ScalarFunction<double>& f) const override;   //!< integral f d3r

    // --- 1E nuclear + external pseudopotential (own the atom/model data). ---
    virtual chmat_t MakeNuclear (const Structure*) const override;  //!< Bare-Coulomb structure factor.

    //! \brief Assemble any local external potential \f$\langle G|V|G'\rangle=\frac1\Omega\sum_a
    //! v(Z_a,|\Delta G|^2)e^{-i\Delta G\cdot\tau_a}\f$, \f$\Delta G\ne 0\f$ (\f$\Delta G=0\f$ dropped).
    virtual chmat_t MakeLocalPotential(const Structure* cl, const Pseudopotential::LocalPotential& loc) const override;
    //! \brief Assemble the separable (Kleinman-Bylander) nonlocal potential (rank-1 per atom, projector, m).
    virtual chmat_t MakeSeparablePotential(const Structure* cl, const Pseudopotential::SeparablePotential& v) const override;

    virtual std::string Name      () const override {return "PlaneWave";}
    virtual std::string BasisSetID() const override; // geometry-aware cache key (Name + k, Ecut, nG)

    virtual std::ostream& Write(std::ostream&) const override;
};

} //namespace
