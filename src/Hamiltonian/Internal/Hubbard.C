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
// Increment 1 carries ONE U per manifold (U_i == U: the shell-averaged special case, exactly what CP2K's
// &DFT_PLUS_U computes); the per-eigenstate / per-site-irrep U_i is the same code with a vector.
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
import qchem.BasisSet.Orbital_DFT_IBS;          // Orbital_DFT_IBS<U,dcmplx> (what ProjectOnto hands the vendor)
import qchem.Mesh.Integrator;                   // qcMesh::MatrixForward / MatrixAdjoint
import qchem.Symmetry;                          // sym_t, SymMap
import qchem.Symmetry.Molecule.ShellRep;        // ShellRep::Rep -- the site group on a shell
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
    //! \a S is the block's overlap; \a columns[M] the block's function indices of manifold M.
    LowdinProjector(const hmat_t<U>& S, const std::vector<std::vector<size_t>>& columns);

    using qcMesh::MatrixForward<U>::Forward;   // un-hide the factored overload (qchem.Mesh.Integrator)
    virtual rvec_t    Forward(const hmat_t<U>& D) const override;
    virtual hmat_t<U> Adjoint(const rvec_t& W)    const override;
    //! The "integral" of an occupation-shaped vector is its TRACE per manifold, summed: the manifold charge.
    virtual double    Integrate(const rvec_t& f) const override;
    virtual size_t    NumCoefficients() const override {return itsN;}
    size_t NumManifolds() const {return itsT.size();}
    size_t Size(size_t M) const {return itsT[M].columns();}   //!< \f$m_M=2l+1\f$
private:
    std::vector<mat_t<U>> itsT;   //!< per manifold: \f$S^{1/2}[:,M]\f$, \f$n\times m\f$
    size_t                itsN;   //!< \f$\sum_M m_M^2\f$
};

//! \brief The DFT+U term.  Periodic (Bloch, dcmplx run) with the real-TRIM corner, spin-native.
class Hubbard_U
    : public virtual cDynamic_HT
    , private        cDynamic_HT_Imp
    , public         Dynamic_HT_RealBlock_Imp
    , public virtual Fitting::ScalarProjector      //!< the FORWARD vendor the density projects onto
{
public:
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
    //! \name Orbital resolution (increment 2): the site group's action on the manifold and its irrep labels
    //!@{
    //! ONE irrep of the manifold's representation, as the term names it: its characters over the site
    //! ops and its dimension.  Two clusters with the same characters are the same irrep.
    struct IrrepSig { rvec_t chi; size_t dim=0; };
    //! ONE LABEL a U can be attached to: (site irrep, dominant grey parent).  Three levels are three labels.
    struct Label { size_t irrep; size_t parent; double parentWeight; size_t dim; };
    struct ManifoldSym
    {
        std::vector<rmat_t>   D;        //!< the site group on the manifold (block-diagonal over its shells)
        std::vector<rmat_t>   Dgrey;    //!< the grey (parent) stabiliser, for parentage only
        std::vector<IrrepSig> irreps;   //!< the site irreps the manifold contains, in a fixed order
        std::vector<IrrepSig> grey;     //!< the grey irreps (t2g, e_g on a cubic site)
        std::vector<rmat_t>   Psite;    //!< the site isotypic projectors (d_k/|G|) sum chi_k(g) D(g)
        std::vector<rmat_t>   Pgrey;    //!< the grey isotypic projectors, for parentage
        std::vector<Label>    labels;   //!< the U slots, printed once at PrepareSlots
        double                purity=1; //!< min over eigenvectors of the chosen irrep's weight (1 = symmetric n)
    };
    //! The manifold's U slots (increment 2): one per (site irrep, grey parent) pair, in \c Uirrep order.
    const std::vector<Label>& Labels(size_t M) const {return itsSym.at(M).labels;}
    //! \f$\sum\lambda\f$ per label of manifold \a M in channel \a s (the occupation of each level).
    const rvec_t& OccupationByLabel(size_t M, const Spin& s) const;

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
    //! The manifold's function indices in \a orb (from its AoShellSource face); throws on a Cartesian d.
    std::vector<std::vector<size_t>> Columns(const BasisSet::Orbital_1E_IBS<double>&) const;
    std::vector<std::vector<size_t>> Columns(const BasisSet::Orbital_1E_IBS<dcmplx>&) const;
    //! The shells' angular reps, one per selected shell in \c Columns order (the symmetry table's input).
    std::vector<std::vector<std::shared_ptr<const Symmetry::Molecule::ShellRep>>> ShellReps(const BasisSet::Orbital_1E_IBS<double>&) const;
    std::vector<std::vector<std::shared_ptr<const Symmetry::Molecule::ShellRep>>> ShellReps(const BasisSet::Orbital_1E_IBS<dcmplx>&) const;
    //! Bring the occupations up to \a cd (a no-op when frozen or already at this serial).
    void EnsureOccupations(const cChargeDensity* cd) const;
    //! Per channel: eigen-decompose \a n (flattened, all manifolds), fill \a occ / \a W, return its \f$E_U\f$.
    double Analyse(const rvec_t& n, std::vector<rvec_t>& occ, rvec_t& W, std::vector<rvec_t>* byLabel=nullptr) const;

    //! \name Orbital resolution -- the private half
    //!@{
    //! Build \c itsSym[M] from the manifold's shells and its \c siteOps (called once, in PrepareSlots).
    void BuildSymmetry(size_t M, const std::vector<std::shared_ptr<const Symmetry::Molecule::ShellRep>>& reps) const;
    //! Eigen-decompose \a n (NEVER symmetrised: the functional is evaluated on the density's own
    //! occupations, symmetry only NAMES them) and label each eigenvector by the site irrep whose isotypic
    //! projector carries most of it, with its grey parent likewise; returns per eigenvalue its label index.
    std::vector<size_t> LabelEigenvectors(size_t M, const rsmat_t& n, rvec_t& lam, rmat_t& v) const;
    const Label& LabelOf(size_t M, size_t label) const {return itsSym[M].labels[label];}
    //!@}

    std::shared_ptr<const Structure> itsSt;
    std::vector<HubbardManifold>     itsManifolds;
    std::vector<rvec3_t>             itsSites;          //!< the cell's site positions, ForEachSite order
    SpinGroup                        itsGroup;
    bool                             itsEigenForm;      //!< RunPolicy::HubbardEigen at construction (false = CP2K's populations)
    mutable size_t                   itsNCoeff = 0;     //!< \f$\sum_M m_M^2\f$, fixed by the first block seen
    bool                             itsFrozen = false;

    struct Channel
    {
        rvec_t              n;     //!< the flattened occupation blocks (the forward's output, k-summed)
        rvec_t              W;     //!< the flattened potential blocks (the adjoint's input)
        std::vector<rvec_t> occ;   //!< per manifold: the eigenvalues \f$\lambda_i\f$
        std::vector<rvec_t> byLabel;   //!< per manifold: \f$\sum\lambda\f$ per label (the printed resolution)
    };
    mutable std::vector<ManifoldSym> itsSym;      //!< per manifold, built in PrepareSlots
    mutable std::map<Spin,Channel> itsChannels;   //!< one per spin irrep of the subgroup
    mutable size_t itsOccVersion = size_t(-1);
    mutable double itsEU = 0.0;
    mutable SymMap<LowdinProjector<dcmplx>> itsProj;    //!< per Bloch block
    mutable SymMap<LowdinProjector<double>> itsProjR;   //!< per real TRIM block
};

} // namespace
