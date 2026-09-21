// File: Hamiltonian/Internal/Hubbard.C  DFT+U -- the Hubbard on-site correction, ORBITAL-RESOLVED by
// construction, born on the MatrixForward/MatrixAdjoint pair.  Programme step 5 (doc/OpenWork.md section 1);
// the ruling is doc/Pins.md pin 23.
//
// WHAT IT IS.  For each Hubbard MANIFOLD M -- the 2l+1 functions of shell l on one atom, an INPUT list on
// the run (never "the transition metal's d shell" by assumption) -- and each spin channel sigma, the
// Loewdin occupation matrix
//     n^sigma_M = Sum_k T_k^dagger D^sigma_k T_k ,      T_k = S_k^{1/2}[:, M]   (n x m, m = 2l+1)
// is eigen-decomposed, n = Sum_i lambda_i v_i v_i^T, and the Dudarev functional in that eigenbasis is
//     E_U = Sum_{sigma,M} Sum_i (U_i/2) lambda_i (1 - lambda_i)                (Macke et al. 2024, eq 6)
//     W^sigma_M = Sum_i U_i (1/2 - lambda_i) v_i v_i^T,   V^sigma_k = T_k W^sigma_M T_k^dagger.
// Increment 1 carried ONE U per manifold (U_i == U: the shell-averaged special case, exactly what CP2K's
// &DFT_PLUS_U computes).  Increment 2 (2026-09-21) is the per-site-irrep U_i: ManifoldSymmetry (below) names
// every eigenvector of the density's OWN n by a (site irrep, chemical parent irrep) SLOT -- a1g<t2g,
// e_g<e_g, e_g<t2g on an AFM-II Mn -- and HubbardManifold::Uirrep is one U per slot.  n is never
// symmetrised; equal Uirrep is bit-identical with the shell-averaged run (gates: UTHamiltonian
// ManifoldSymmetry.*, UTStructure SiteGroups.*, gpwprobe mno MNO_U_IRREP).
//
// ★ BORN ON THE PAIR (user, 2026-09-11).  The occupation-matrix FORWARD (D -> n) and the potential ADJOINT
// (W -> V) are the two faces of ONE object per orbital block, LowdinProjector<TBlock>: the same T_k
// contracted both ways, so H = dE/dD holds exactly.  The density gets the forward through its own
// ProjectOnto (this term realises Fitting::ScalarProjector for it, exactly as DeltaScalarFitter does on
// the XC route); the term keeps the adjoint.  Neither can name the other direction.
//
// ★ SCALAR-GENERIC (ruled 2026-09-19, so CleanupCandidates V1.35 collapses it for free): ONE template body
// per operation over the BLOCK scalar; MakeMatrix / MakeMatrixR are one-line instantiations, as Vee_Hartree.
//
// WHERE D COMES FROM.  The Fock build sees the MIXED density, which on the Kerker/Pulay recipes has no D;
// the occupations are taken from the DM-backed SOURCE the mixer retains (tDM_Sourced_CD::DMSource, one cast
// away, as XC's cusp-deficit route does).  A matrix-free SEED gives n = 0 (V = U/2 P, E_U = 0) -- CP2K
// ramps U from zero for the same reason.  Declared: "+U occupations from D_out" is a CP2K deviation (they
// mix P); same fixed point, different trajectory.  SpinGroup::None: n = n_tot/2 in both channels -- the
// zeta=0 collapse (pin 5).
//
// THE FORM (RunPolicy::HubbardEigen, knob QCHEM_U_EIGEN, off under CP2K_COMPAT): Dudarev on the block's
// EIGENVALUES (above) or CP2K's DIAGONAL POPULATIONS -- dft_plus_u.F keeps only q_ii, so E_U = U/2 Sum q_ii(1-q_ii)
// and V_ii = U(1/2 - q_ii), no eigen-decomposition.  Read once at construction, like every deviation.
//
// FROZEN-OCCUPATION MODE (day one): FreezeOccupations(true) keeps the stored n^sigma through every refresh
// -- Macke fixes V_U at its unperturbed value during LR-cDFT perturbation runs, and a polaron study puts the
// carrier on a chosen site the same way.
module;
#include <complex>
#include <iosfwd>
#include <map>
#include <memory>
#include <string>
#include <vector>
export module qchem.Hamiltonian.Internal.Hubbard;
import qchem.Hamiltonian.Internal.Term;        // cDynamic_HT + the _Imp cache mixins (Bloch + real TRIM)
export import qchem.Hamiltonian.Factory;        // HubbardManifold (the public input vocabulary)
import qchem.Hamiltonian.Types;                 // cobs_t / robs_t / tobs_t<U>
import qchem.Fitting.FunctionFitter;            // Fitting::ScalarProjector (the forward vendor face)
import qchem.BasisSet.Orbital_DFT_IBS;
import qchem.BasisSet.BareCoulombSource;         // ERI4Block (ManifoldIntegrals)          // Orbital_DFT_IBS<U,dcmplx> (what ProjectOnto hands the vendor)
import qchem.Mesh.Integrator;                   // qcMesh::MatrixForward / MatrixAdjoint
import qchem.Symmetry;                          // sym_t, SymMap
import qchem.Symmetry.Molecule.OperationRep;    // AoShell (rep + per-component norm) -- the site group on a shell
import qchem.Symmetry.Irrep;                    // Spin, SpinGroup, SpinIrreps
import qchem.Structure;
import qchem.Blaze;

export namespace qchem::Hamiltonian
{

//! \brief THE PAIR for one orbital block: \f$n=T^\dagger DT\f$ forward, \f$V=TWT^\dagger\f$ adjoint,
//! \f$T=S^{1/2}[:,M]\f$ per manifold.  \c NumCoefficients is \f$\sum_M m_M^2\f$: the real
//! \f$m\times m\f$ occupation blocks of every manifold, flattened row-major and concatenated -- the
//! coefficient axis the term reads and writes.
template <class U> class LowdinProjector
    : public virtual qcMesh::MatrixForward<U>
    , public virtual qcMesh::MatrixAdjoint<U>
{
public:
    //! \a S is the block's overlap; \a columns[M] the block's function indices of manifold M; \a contraction[M]
    //! is EMPTY (the manifold's functions ARE those columns: \f$T=S^{1/2}[:,M]\f$, CP2K's convention) or a
    //! \f$|{\rm cols}_M|\times m_M\f$ matrix \f$V\f$ whose columns are the manifold's functions as combinations
    //! of those columns (a CONTRACTED radial, increment 3) -- the pair then S-orthonormalises them within the
    //! manifold, \f$\tilde V=V\,(V^\dagger S_{cc}V)^{-1/2}\f$, and projects onto THOSE functions: the ATOMIC
    //! projector \f$T=S[:,{\rm cols}]\tilde V\f$ (\f$T^\dagger c=\langle\chi|\psi\rangle\f$, QE's `atomic`), where a
    //! column manifold is LÖWDIN (\f$S^{1/2}\f$) -- see the ctor for why the two are not one formula.
    //! \a ortho[M] (contracted manifolds only): join the ORTHO-ATOMIC set -- after each manifold's own
    //! S-orthonormalisation, the functions of every flagged manifold are Löwdin-orthogonalised among
    //! themselves, \f$\tilde W=W\,(W^\dagger SW)^{-1/2}\f$, and projected: \f$T_M=S\tilde W_M\f$.  The
    //! bare-integral map \c Contraction(M) keeps the on-site \f$\tilde V_M\f$.
    LowdinProjector(const hmat_t<U>& S, const std::vector<std::vector<size_t>>& columns,
                    const std::vector<mat_t<U>>& contraction = {}, const std::vector<bool>& ortho = {});

    using qcMesh::MatrixForward<U>::Forward;   // un-hide the factored overload (qchem.Mesh.Integrator)
    virtual rvec_t    Forward(const hmat_t<U>& D) const override;
    virtual hmat_t<U> Adjoint(const rvec_t& W)    const override;
    //! The "integral" of an occupation-shaped vector is its TRACE per manifold, summed: the manifold charge.
    virtual double    Integrate(const rvec_t& f) const override;
    virtual size_t    NumCoefficients() const override {return itsN;}
    size_t NumManifolds() const {return itsT.size();}
    size_t Size(size_t M) const {return itsT[M].columns();}   //!< \f$m_M=2l+1\f$
    const mat_t<U>& T(size_t M) const {return itsT[M];}       //!< the manifold's \f$S^{1/2}\tilde V\f$ (the ACBN0 estimator's Löwdin coefficients \f$\ell=T^\dagger c\f$)
    //! THE COEFFICIENT MAP \f$Q_M\f$ (\f$m\times n\f$): an orbital's AO coefficients \f$c\f$ → its coefficients on the
    //! manifold's functions, \f$Q_Mc\f$.  Column manifold: the selector (the raw coefficients on those columns --
    //! the ACBN0 paper's eq 9); contracted manifold: \f$\tilde V^\dagger S[{\rm cols},:]\f$, the S-metric
    //! projection onto the S-orthonormal \f$\chi_m\f$ (exact for a function in their span).
    const mat_t<U>& Coefficients(size_t M) const {return itsQ[M];}
    //! The manifold's functions over its columns (\f$\tilde V\f$; the identity selector for a column manifold) --
    //! what carries the columns' bare integrals to the manifold's (ERI4Block::Transform).
    const rmat_t& Contraction(size_t M) const {return itsV[M];}
    const std::vector<size_t>& Columns(size_t M) const {return itsCols[M];}
private:
    std::vector<mat_t<U>> itsT;   //!< per manifold: \f$S^{1/2}\tilde V\f$, \f$n\times m\f$
    std::vector<mat_t<U>> itsQ;   //!< per manifold: the coefficient map, \f$m\times n\f$
    std::vector<rmat_t>   itsV;   //!< per manifold: \f$\tilde V\f$ over the columns (real: a radial contraction)
    std::vector<std::vector<size_t>> itsCols;
    size_t                itsN;   //!< \f$\sum_M m_M^2\f$
};

//! \brief THE SITE GROUP ON ONE HUBBARD MANIFOLD, and the names it gives the occupation eigenvectors
//! (increment 2, the orbital resolution of pin 23).
//!
//! The manifold's \f$m\f$ functions carry a representation \f$D(g)\f$ of the site group (block-diagonal
//! over its shells, one \c ShellRep::Rep block each).  Its irreps are CLUSTERED, not looked up -- the
//! molecular character tables are abelian-only and a d shell on a cubic or trigonal site carries 3-D and
//! 2-D irreps: a generic symmetric matrix symmetrised under \f$D\f$ has degenerate clusters that ARE irrep
//! copies, and a cluster's character vector \f$\chi(g)=\mathrm{Tr}(P\,D(g))\f$ names it.  The isotypic
//! projectors \f$P_k\propto\sum_g\chi_k(g)D(g)\f$ then say, for ANY vector, how much of it lies in irrep k.
//!
//! TWO GROUPS.  The SITE group is the declared magnetic decoration's stabiliser (AFM-II MnO: \f$D_{3d}\f$,
//! which splits \f$t_{2g}\to a_{1g}+e_g\f$); the GREY group is the spin-blind stabiliser (\f$O_h\f$), kept
//! for PARENTAGE only -- so a level is named "\f$e_g\f$ (site) < \f$t_{2g}\f$ (grey)".  A U SLOT is the
//! INTERSECTION of a site isotypic component with a grey one, \f$\dim=\mathrm{Tr}(P^{site}_kP^{grey}_p)\f$
//! (the two commute: the site group is a subgroup), fixed by group theory before any density -- so
//! \c HubbardManifold::Uirrep addresses slots that never move.  MnO's d under \f$D_{3d}<O_h\f$: {1, 2, 2}.
//!
//! \f$n\f$ IS NEVER SYMMETRISED (user, 2026-09-21): the functional is evaluated on the density's OWN
//! occupations; symmetry only NAMES their eigenvectors -- each takes the slot whose projectors carry most
//! of it, and the labelling reports its \c purity (min site weight) and per-slot \c parentage (min grey
//! weight) so an impure n is visible, never hidden.  Inside a DEGENERATE cluster the eigenbasis is
//! LAPACK's arbitrary choice, so the cluster is rotated to diagonalise the projectors first (the seed
//! n = 0 is one big cluster; an exact symmetric n has them by construction) -- the spectral decomposition
//! of n is unchanged, the names become well-defined.
class ManifoldSymmetry
{
public:
    //! ONE irrep of the manifold's representation: its characters over the ops and its dimension.  Two
    //! clusters with the same characters are the same irrep.
    struct IrrepSig { rvec_t chi; size_t dim=0; };
    //! ONE U SLOT: (site irrep, grey parent) with \f$\dim=\mathrm{Tr}(P^{site}_kP^{grey}_p)\f$ -- more than
    //! the irrep's dimension when the irrep occurs several times inside one parent.
    struct Slot { size_t irrep; size_t parent; size_t dim; };
    //! What one eigen-decomposition was told about itself.
    struct Labelling
    {
        std::vector<size_t> slot;        //!< per eigenvalue: its slot index
        double              purity=1.0;  //!< min over eigenvectors of the chosen site irrep's weight (1 = symmetric n)
        rvec_t              parentage;   //!< per slot: min over its eigenvectors of the grey parent's weight (1 if unused)
    };

    //! \a D the site group on the manifold, \a Dgrey the grey stabiliser on it (the same list when the run
    //! has no magnetic decoration).  Both lists must contain the identity and be closed under products.
    ManifoldSymmetry(std::vector<rmat_t> D, std::vector<rmat_t> Dgrey);
    //! The manifold representation of the Cartesian ops \a ops on the shells \a shells: block-diagonal, one
    //! block per shell in \a shells order, each \f$(N_a/N_b)\,\mathrm{Rep}(b,a)\f$ -- the shell's angular
    //! rep carried from its RAW components to its NORMALISED ones by \c AoShell::norm, exactly as
    //! \c BuildOperationRep does (the occupation matrix lives in the normalised functions; the raw real
    //! harmonics' norms differ across m, so the raw rep is not orthogonal).  An empty op list is
    //! \f$C_1\f$ (the identity only).
    static std::vector<rmat_t> Rep(const std::vector<Symmetry::Molecule::AoShell>& shells, const std::vector<rmat3d_t>& ops);

    //! Eigen-decompose \a n (\a lam ascending, \a v its eigenvectors in columns) and name each eigenvector.
    Labelling Label(const rsmat_t& n, rvec_t& lam, rmat_t& v) const;

    size_t Size()       const {return itsM;}              //!< \f$m\f$, the manifold's function count
    size_t NumOps()     const {return itsD.size();}
    size_t NumGreyOps() const {return itsDgrey.size();}
    bool   Graded()     const {return itsDgrey.size()!=itsD.size();}   //!< is there a parentage to report?
    const std::vector<IrrepSig>& Irreps() const {return itsIrreps;}
    const std::vector<IrrepSig>& Grey()   const {return itsGrey;}
    const std::vector<Slot>&     Slots()  const {return itsSlots;}
    const rmat_t& SiteProjector(size_t k) const {return itsPsite[k];}
    const rmat_t& GreyProjector(size_t p) const {return itsPgrey[p];}
    //! The slot table as the term prints it: "[k] dim d (site irrep i) < grey irrep p dim dp  U=... eV".
    std::ostream& WriteSlots(std::ostream&, double U, const std::vector<double>& Uirrep) const;

private:
    size_t                itsM;
    std::vector<rmat_t>   itsD, itsDgrey;
    std::vector<IrrepSig> itsIrreps, itsGrey;
    std::vector<rmat_t>   itsPsite, itsPgrey;
    std::vector<Slot>     itsSlots;
};

//! \brief Dudarev in the occupation eigenbasis with ONE U PER SLOT (Macke et al. 2024 eq 6):
//! \f$E_U=\sum_i\tfrac{U_i}2\lambda_i(1-\lambda_i)\f$, \f$W=\sum_iU_i(\tfrac12-\lambda_i)v_iv_i^T\f$
//! rotated back to the manifold's own basis.  \f$U_i\f$ is \a Uirrep[slot[i]] when \a Uirrep is filled
//! (and \a slot names the eigenvector), else \a U -- the shell-averaged special case.  Returns \f$E_U\f$.
double DudarevInEigenbasis(const rvec_t& lam, const rmat_t& v, const std::vector<size_t>& slot,
                           const std::vector<double>& Uirrep, double U, rmat_t& W);

//! \brief WHAT A U-ESTIMATOR CONSUMES FROM THE +U TERM (increment 3, ACBN0): the manifolds, their
//! equivalence (same species and l -- the paper's \f$\{\bar m\}\f$ set the renormalised occupation is
//! summed over), each block's manifold functions (for the bare integrals) and the Löwdin coefficients of
//! given orbitals in every manifold.  An abstract face so the composite Hamiltonian finds the term by an
//! abstract->abstract cast and the estimator never names \c Hubbard_U.
class HubbardProjection
{
public:
    virtual ~HubbardProjection() = default;
    virtual const std::vector<HubbardManifold>& Manifolds() const = 0;
    virtual SpinGroup Group() const = 0;
    //! The manifolds of the same (species, l) as \a M, \a M itself included.
    virtual std::vector<size_t> EquivalentManifolds(size_t M) const = 0;
    //! Per manifold: the BARE two-electron integrals over the manifold's functions on \a block (a column
    //! manifold: the block's own; a contracted one: carried through the contraction on all four indices).
    virtual std::vector<BasisSet::ERI4Block> ManifoldIntegrals(const BasisSet::Orbital_DFT_IBS<double,dcmplx>& block) const = 0;
    virtual std::vector<BasisSet::ERI4Block> ManifoldIntegrals(const BasisSet::Orbital_DFT_IBS<dcmplx,dcmplx>& block) const = 0;
    //! Per manifold: the orbitals' COEFFICIENTS on the manifold's functions, \f$Q_MC\f$ (\f$m_M\times n_{\rm orb}\f$;
    //! \c LowdinProjector::Coefficients) -- what pairs with \c ManifoldIntegrals in an on-site HF energy.
    virtual std::vector<mat_t<double>> ManifoldCoefficients(const BasisSet::Orbital_DFT_IBS<double,dcmplx>& block, const mat_t<double>& C) const = 0;
    virtual std::vector<mat_t<dcmplx>> ManifoldCoefficients(const BasisSet::Orbital_DFT_IBS<dcmplx,dcmplx>& block, const mat_t<dcmplx>& C) const = 0;
    //! Per manifold: \f$\ell=T_M^\dagger C\f$ (\f$m_M\times n_{\rm orb}\f$), the LÖWDIN coefficients -- what a charge is made of.
    virtual std::vector<mat_t<double>> LowdinCoefficients(const BasisSet::Orbital_DFT_IBS<double,dcmplx>& block, const mat_t<double>& C) const = 0;
    virtual std::vector<mat_t<dcmplx>> LowdinCoefficients(const BasisSet::Orbital_DFT_IBS<dcmplx,dcmplx>& block, const mat_t<dcmplx>& C) const = 0;
    //! THE OTHER DIRECTION -- what the estimator does to the term: set manifold \a M's \f$U\f$ (Hartree; a
    //! filled \c Uirrep is set to it throughout, shell-averaged) for the NEXT Fock build.  The outer loop of
    //! the paper (SCF at \f$U^{(n)}\f$, estimate, run again) lives on this.
    virtual void SetU(size_t M, double U) = 0;
};

//! \brief The DFT+U term.  Periodic (Bloch, dcmplx run) with the real-TRIM corner, spin-native.
class Hubbard_U
    : public virtual cDynamic_HT
    , private        cDynamic_HT_Imp
    , public         Dynamic_HT_RealBlock_Imp
    , public virtual Fitting::ScalarProjector      //!< the FORWARD vendor the density projects onto
    , public virtual HubbardProjection             //!< what the ACBN0 estimator consumes (increment 3)
{
public:
    //! What a manifold is on a block: its column indices (from the AoShellSource face; throws on a Cartesian
    //! d), the raw contraction over them (empty = the columns themselves), and the shells that carry its
    //! angular rep (ALL the selected shells for a column manifold; ONE for a contracted one -- they share it).
    struct Selection
    {
        std::vector<std::vector<size_t>>                      cols;
        std::vector<rmat_t>                                   contraction;
        std::vector<std::vector<Symmetry::Molecule::AoShell>> shells;
    };
    Selection Select(const BasisSet::Orbital_1E_IBS<double>&) const;
    Selection Select(const BasisSet::Orbital_1E_IBS<dcmplx>&) const;

    //! \name HubbardProjection
    //!@{
    virtual const std::vector<HubbardManifold>& Manifolds() const override {return itsManifolds;}
    virtual SpinGroup Group() const override {return itsGroup;}
    virtual std::vector<size_t> EquivalentManifolds(size_t M) const override;
    virtual std::vector<BasisSet::ERI4Block> ManifoldIntegrals(const BasisSet::Orbital_DFT_IBS<double,dcmplx>&) const override;
    virtual std::vector<BasisSet::ERI4Block> ManifoldIntegrals(const BasisSet::Orbital_DFT_IBS<dcmplx,dcmplx>&) const override;
    virtual std::vector<mat_t<double>> ManifoldCoefficients(const BasisSet::Orbital_DFT_IBS<double,dcmplx>&, const mat_t<double>&) const override;
    virtual std::vector<mat_t<dcmplx>> ManifoldCoefficients(const BasisSet::Orbital_DFT_IBS<dcmplx,dcmplx>&, const mat_t<dcmplx>&) const override;
    virtual std::vector<mat_t<double>> LowdinCoefficients(const BasisSet::Orbital_DFT_IBS<double,dcmplx>&, const mat_t<double>&) const override;
    virtual std::vector<mat_t<dcmplx>> LowdinCoefficients(const BasisSet::Orbital_DFT_IBS<dcmplx,dcmplx>&, const mat_t<dcmplx>&) const override;
    virtual void SetU(size_t M, double U) override;
    //!@}

    //! \a st names the sites the manifolds index; \a g the imposed spin subgroup (which channels exist).
    Hubbard_U(const std::shared_ptr<const Structure>& st, std::vector<HubbardManifold> manifolds, SpinGroup g);
    ~Hubbard_U();

    //! \copydoc HT_SlotOwner::PrepareSlots  (two caches: the Bloch one and the real TRIM one) -- AND the
    //! geometry phase of this term: every block's pair is built here, so \c NumCoefficients is known before
    //! the first density asks for its occupations (the composite sizes its sum by it BEFORE any block runs).
    virtual void PrepareSlots(const cbs_t* bs) const override;
    //! THE EAGER PHASE: the occupations are k-independent, so they are taken here, once per density.
    virtual void          RefreshForDensity(const cChargeDensity* cd) const override;
    virtual void          GetEnergy(EnergyBreakdown&, const cDM_CD*) const override;
    virtual std::ostream& Write(std::ostream&) const override;
    virtual bool          IsVirialValid() const override {return false;}   //!< an occupation functional is not Coulombic

    //! Hold the current occupations through every refresh (LR-cDFT perturbation runs; occupation control).
    void FreezeOccupations(bool frozen) {itsFrozen=frozen;}
    bool OccupationsFrozen() const {return itsFrozen;}
    //! The current occupation eigenvalues of manifold \a M, channel \a s (empty before the first refresh).
    const rvec_t& Occupations(size_t M, const Spin& s) const;
    double HubbardEnergy() const {return itsEU;}   //!< \f$E_U\f$ of the last refresh
    //! \name Orbital resolution (increment 2): the U slots and what each level holds
    //!@{
    //! The manifold's U slots: one per (site irrep, grey parent) intersection, in \c Uirrep order; empty
    //! before \c PrepareSlots (or under CP2K's population form, which has no eigenvectors to name).
    const std::vector<ManifoldSymmetry::Slot>& Slots(size_t M) const;
    //! \f$\sum\lambda\f$ per slot of manifold \a M in channel \a s (the occupation of each level).
    const rvec_t& OccupationByLabel(size_t M, const Spin& s) const;
    //! How well the last density's eigenvectors sat in the slots (1 = a site-symmetric \f$n\f$).
    const ManifoldSymmetry::Labelling& LastLabelling(size_t M, const Spin& s) const;
    //!@}

    //! \name Fitting::ScalarProjector -- the forward vendor (what tDM_CD::ProjectOnto asks)
    //!@{
    virtual const qcMesh::MatrixForward<double>& Forward(const BasisSet::Orbital_DFT_IBS<double,dcmplx>&) const override;
    virtual const qcMesh::MatrixForward<dcmplx>& Forward(const BasisSet::Orbital_DFT_IBS<dcmplx,dcmplx>&) const override;
    //! A matrix-free density has no occupation matrix: zeros (a seed; see the header).
    virtual rvec_t Project(const ScalarFunction<double>&) const override {return rvec_t(NumCoefficients(),0.0);}
    virtual size_t NumCoefficients() const override {return itsNCoeff;}
    //!@}

private:
    virtual chmat_t MakeMatrix (const cobs_t*, const Spin&, const cChargeDensity*) const override;
    virtual rsmat_t MakeMatrixR(const robs_t*, const Spin&, const cChargeDensity*) const override;
    template <class U> hmat_t<U> MakeMatrixT(const tobs_t<U>*, const Spin&, const cChargeDensity*) const;

    //! The block's pair, built on first use (the overlap is geometry-fixed) and keyed like the XC route's.
    template <class U> const LowdinProjector<U>& Projector(const BasisSet::Orbital_DFT_IBS<U,dcmplx>&) const;
    template <class U> SymMap<LowdinProjector<U>>& Cache() const;
    //! Bring the occupations up to \a cd (a no-op when frozen or already at this serial).
    void EnsureOccupations(const cChargeDensity* cd) const;
    //! Per channel: eigen-decompose \a n (flattened, all manifolds), fill the channel's occupations, slot
    //! sums, labelling and \f$W\f$, return its \f$E_U\f$.
    struct Channel;
    double Analyse(Channel& ch) const;

    //! Build \c itsSym[M] from the manifold's shells and its \c siteOps / \c greyOps (once, in PrepareSlots),
    //! and say the slot table once (pin 17).
    void BuildSymmetry(size_t M, const std::vector<Symmetry::Molecule::AoShell>& shells) const;

    std::shared_ptr<const Structure> itsSt;
    std::vector<HubbardManifold>     itsManifolds;
    std::vector<rvec3_t>             itsSites;          //!< the cell's site positions, ForEachSite order
    std::vector<int>                 itsSiteZ;          //!< their species (EquivalentManifolds)
    SpinGroup                        itsGroup;
    bool                             itsEigenForm;      //!< RunPolicy::HubbardEigen at construction (false = CP2K's populations)
    mutable size_t                   itsNCoeff = 0;     //!< \f$\sum_M m_M^2\f$, fixed by the first block seen
    bool                             itsFrozen = false;

    struct Channel
    {
        rvec_t              n;     //!< the flattened occupation blocks (the forward's output, k-summed)
        rvec_t              W;     //!< the flattened potential blocks (the adjoint's input)
        std::vector<rvec_t> occ;   //!< per manifold: the eigenvalues \f$\lambda_i\f$
        std::vector<rvec_t> byLabel;   //!< per manifold: \f$\sum\lambda\f$ per slot (the printed resolution)
        std::vector<ManifoldSymmetry::Labelling> lab;   //!< per manifold: the last labelling's diagnostics
    };
    mutable std::vector<std::shared_ptr<const ManifoldSymmetry>> itsSym;   //!< per manifold, built in PrepareSlots (null = C_1 or the population form)
    mutable std::map<Spin,Channel> itsChannels;   //!< one per spin irrep of the subgroup
    mutable size_t itsOccVersion = size_t(-1);
    mutable double itsEU = 0.0;
    mutable SymMap<LowdinProjector<dcmplx>> itsProj;    //!< per Bloch block
    mutable SymMap<LowdinProjector<double>> itsProjR;   //!< per real TRIM block
};

} // namespace
