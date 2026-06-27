// File: BasisSet/Lattice_3D/PlaneWave_IBS.C  Plane-wave irrep basis set for one k-point.
//
// A complex (dcmplx) Orbital_1E_IBS whose functions are the normalised plane waves
// e^{i(k+G).r}/sqrt(V) for the reciprocal lattice vectors G in the cutoff set
// { G : 1/2 |k+G|^2 < Ecut }.  The wave-vector k labels the Bloch (translational)
// symmetry of the block (BlochFactory).
//
// Milestone 1 (empty lattice, V_ext=0): only Overlap (identity) and Kinetic
// (diagonal |k+G|^2) are exercised; the eigenvalues of 1/2*Kinetic reproduce the
// exact free-electron ladder 1/2 |k+G|^2.  See doc/PlaneWavePlan.md (sections 2.1, 3, 6).
module;
#include <functional>
#include <iosfwd>
#include <string>
#include <vector>

export module qchem.BasisSet.Lattice_3D.PlaneWave_IBS;
export import qchem.BasisSet.Band_DFT_IBS;     // the abstract real-space DFT-integration capability
export import qchem.BasisSet.Band_FT_IBS;       // the abstract G-space DFT capability (+ FourierMap)
import qchem.BasisSet.Internal.IrrepBasisSetImp;   // IrrepBasisSetImp<T>: GetSymmetry/GetSymt/GetIrrep
export import qchem.ReciprocalLattice;             // ctor takes a ReciprocalLattice (carries the B cell)
export import qchem.Pseudopotential.LocalPotential;      // local potential form-factor abstraction
export import qchem.Pseudopotential.SeparablePotential; // KB nonlocal projector abstraction
import qchem.Structure;
import qchem.Symmetry;                             // sym_t (the Bloch irrep handed to the ctor)
import qchem.Types;

export namespace BasisSet::Lattice_3D
{

//! \brief Plane-wave basis for a single k-point: the normalised waves
//! \f$ e^{i(k+G)\cdot r}/\sqrt V \f$ over the cutoff set \f$\{G:\tfrac12|k+G|^2<E_{cut}\}\f$.
class PlaneWave_IBS
    : public virtual BasisSet::Band_DFT_IBS<dcmplx> // real-space DFT-integration (Hartree/XC/external)
    , public virtual BasisSet::Band_FT_IBS           // G-space DFT (rho-tilde -> Hartree, FFT XC)
    , public         BasisSet::IrrepBasisSetImp<dcmplx> // supplies GetSymmetry/GetSymt/GetIrrep + itsSymmetry
{
public:
    //! \brief Primary constructor: the Bloch symmetry IS the k-label (mirrors the atom IBSs, which take
    //! an abstract \c sym_t and pry out their quantum number).  The crystal momentum is read from the
    //! irrep via Symmetry::Getk; the basis owns no copy of the BZ grid.
    //! \param recip   the reciprocal lattice (its UnitCell matrix is \f$B=2\pi A^{-\top}\f$).
    //! \param irrep   the Bloch irrep (a BlochQN); \f$k\f$ = Symmetry::Getk(irrep).
    //! \param Ecut    plane-wave energy cutoff (Hartree): keep \f$G\f$ with \f$\tfrac12|k+G|^2<E_{cut}\f$.
    PlaneWave_IBS(const ReciprocalLattice& recip, const sym_t& irrep, double Ecut);

    //! \brief Convenience constructor for tests/callers that work directly in BZ-grid indices: builds
    //! the Bloch irrep \c BlochFactory(N,kIndex) and delegates to the primary constructor above.
    //! \param N       Brillouin-zone grid divisions (context for the integer k-label).
    //! \param kIndex  integer k-label; the fractional crystal momentum is \f$k = kIndex/N\f$.
    PlaneWave_IBS(const ReciprocalLattice& recip, const ivec3_t& N,
                  const ivec3_t& kIndex, double Ecut);

    virtual size_t GetNumFunctions() const {return itsG.size();}

    //! Reciprocal-index label \f$m\f$ of basis function \a i (its plane wave is \f$e^{i(k+G)\cdot r}\f$,
    //! \f$G = B\,m\f$).  This integer triple is the defining quantum number of the plane wave.
    virtual ivec3_t GetGIndex(size_t i) const {return itsG[i];}

    // --- Band_DFT_IBS capability: high-level "do XYZ" questions the KS potential terms ask. ---
    // The basis owns the integration (its own FFT grid); the terms never see G-vectors or the mesh.
    virtual chmat_t Overlap  (const ScalarFunction<double>& f) const override;            //!< <i|f|j> (uncached)
    virtual chmat_t Repulsion(const ScalarFunction<double>& rho, double& Eh) const override; //!< <i|V_Coul[rho]|j> (uncached)
    virtual double  Integral         (const ScalarFunction<double>& f) const;            //!< integral f d3r

    //! \brief Density Fourier coefficients \f$\tilde\rho(\Delta m)=\frac1\Omega\sum_{G_i-G_j=\Delta m}D_{ij}\f$
    //! for a density matrix \a D in THIS plane-wave block.  The G-space route to the Hartree/XC matrices:
    //! no \f$O(N_{pts}n^2)\f$ real-space sampling -- one \f$O(n^2)\f$ accumulation over the difference set.
    virtual FourierMap MakeFourierDensity(const chmat_t& D) const override;
    //! \brief Hartree matrix + energy directly from the density's G-space coefficients \a rho
    //! (= MakeFourierDensity): \f$V_H(\Delta m)=4\pi\tilde\rho/|B\Delta m|^2\f$, \f$E_H=\tfrac\Omega2\sum
    //! 4\pi|\tilde\rho|^2/G^2\f$ (\f$\Delta m=0\f$ dropped).  The FFT-free Poisson solve.
    virtual chmat_t Repulsion(const FourierMap& rho, double& Eh) const override;

    // XC route (basis owns the FFTs; see Band_FT_IBS).  rho(r) via inverse FFT of rho-tilde; the
    // term maps the functional over the grid values and hands them back for the forward FFT / quadrature.
    virtual rvec_t     RhoOnGrid   (const FourierMap& rho) const override;
    virtual FourierMap ForwardGrid (const rvec_t& gridValues) const override; //!< forward FFT -> Vtilde(dm)
    virtual chmat_t    Overlap     (const FourierMap& Vtilde) const override;  //!< <i|V|j>=Vtilde(dm)
    virtual chmat_t    Overlap     (const rvec_t& Vgrid)  const override;      //!< = Overlap(ForwardGrid(Vgrid))
    virtual double     Integral    (const rvec_t& fgrid)   const override;
    //!< (N/Omega) Sum_a alpha_a for the supplied local model (the dropped-G=0 alignment energy).
    virtual double  ExternalG0Energy(const Structure* cl, const std::function<double(int)>& formFactorG0, double numElectrons) const override;

    // 1E integral building blocks (no 1/2 on Kinetic -- the Hamiltonian applies it).
    virtual chmat_t MakeOverlap () const;                  //!< Identity (PWs orthonormal over the cell).
    virtual chmat_t MakeKinetic () const;                  //!< Diagonal \f$|k+G|^2 = \langle p^2\rangle\f$.
    virtual chmat_t MakeNuclear (const Structure*) const;  //!< Bare-Coulomb structure factor (= MakeLocalPotential with BareCoulomb).

    //! \brief Assemble any local external potential: \f$ \langle G|V|G'\rangle = \frac1\Omega \sum_a
    //! v(Z_a,|\Delta G|^2)\, e^{-i\Delta G\cdot\tau_a} \f$, \f$\Delta G=G-G'\ne 0\f$ (the \f$\Delta G=0\f$
    //! term is dropped -- uniform neutralising background).  The species form factor \a v selects the
    //! potential model (bare Coulomb, Gaussian-smeared nucleus, pseudopotential, ...).  Hermitian.
    virtual chmat_t MakeLocalPotential(const Structure* cl, const std::function<double(int,double)>& formFactor) const override;

    //! \brief Assemble the separable (Kleinman-Bylander) NONLOCAL potential \f$ \langle G|V_{NL}|G'\rangle
    //! = \frac1\Omega \sum_a e^{-i\Delta G\cdot\tau_a} \sum_p \tilde\beta_p(|k+G|)\,D_p\,\tilde\beta_p(|k+G'|)\f$.
    //! The external one-body potential is then \f$V = V_{loc} + V_{NL}\f$ (both summed into \f$H(k)\f$).
    //! Hermitian; per atom & projector this is a rank-1 \f$|\beta\rangle D\langle\beta|\f$ contribution.
    virtual chmat_t MakeSeparablePotential(const Structure* cl, const Pseudopotential::SeparablePotential& v) const override;

    //! \brief Assemble \f$ \langle G|V|G'\rangle = \tilde V(G-G') \f$ from a caller-supplied G-space
    //! potential, keyed by the reciprocal-index difference \f$\Delta m = m(G)-m(G')\f$.
    //!
    //! This is the reusable G-space potential assembly: the milestone-2.2 separable cosine and the
    //! milestone-2.3 nuclear structure factor are both just particular \f$\tilde V\f$ suppliers.
    //! \note Returns smat_t (symmetric) -- correct for a real, even \f$\tilde V\f$ (the cosine).  A
    //! \f$\tilde V\f$ should satisfy \f$\tilde V(-\Delta m)=\overline{\tilde V(\Delta m)}\f$ so the
    //! result is Hermitian (chmat_t); this holds for the cosine and the nuclear structure factor.
    virtual chmat_t MakePotential(const std::function<dcmplx(const ivec3_t&)>& Vtilde) const;

    // VectorFunction: the plane-wave values / gradients at a point.
    virtual cvec_t     operator() (const rvec3_t& r) const;
    virtual cvec3vec_t Gradient   (const rvec3_t& r) const;

    virtual std::string Name      () const {return "PlaneWave";}
    virtual std::string BasisSetID() const; // geometry-aware cache key (k, Ecut, nG)

    virtual std::ostream& Write(std::ostream&) const;

private:
    rvec3_t GetGCartesian(const ivec3_t& m) const; //!< \f$ G = B\,m \f$ in Cartesian a.u.
    std::vector<rvec3_t> UniformGrid(const ivec3_t& npts) const; //!< fractional XC grid (Omega/prod weight).
    ivec3_t AutoGrid() const;       //!< grid divisions resolving the basis difference set (no aliasing).
    ivec3_t FFTGrid()  const;       //!< AutoGrid padded to powers of two (radix-2 FFT for the XC route).

    ReciprocalLattice    itsRecip;  //!< Reciprocal cell (matrix \f$B\f$); source of G and |k+G|.
    rvec3_t              itsk;      //!< Fractional crystal momentum \f$k\f$ (read from the Bloch irrep).
    double               itsEcut;   //!< Energy cutoff (Hartree).
    double               itsVolume; //!< Direct cell volume \f$V\f$ (for the \f$1/\sqrt V\f$ norm).
    std::vector<ivec3_t> itsG;      //!< Surviving reciprocal index triples \f$m\f$ (\f$G=Bm\f$).
};

} //namespace
