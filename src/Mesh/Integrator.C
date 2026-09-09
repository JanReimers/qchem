// File: Mesh/Integrator.C  MatrixIntegrator -- a FORWARD/ADJOINT PAIR that cannot be mismatched.
//
// ★ WHY A CLASS AND NOT TWO FREE FUNCTIONS (user, 2026-09-08).
//
// A density-functional assembly runs one map in each direction:
//
//   FORWARD  D  ->  rho(r_g)        collocate the density matrix onto the points
//   ADJOINT  v  ->  <i|v|j>         integrate a field on those points back to a matrix
//
// \f$H=\partial E/\partial D\f$ holds only if the adjoint is the EXACT adjoint of the forward ON THE SAME
// TRUNCATED OPERATOR.  Offered as two independent free functions, a mismatch is EXPRESSIBLE -- and this
// tree has measured what that costs: an unscreened \f$\rho\f$ paired with a screened \f$H\f$ sent Si from
// 14 to 60 SCF iterations and moved E by 35 µHa (doc/CleanupCandidates.md R1.0j).  Behind one interface
// the pairing is a class invariant: a caller cannot obtain half of one route and half of another.
//
// ⚠ THE ASSUMPTION THAT MAKES THAT AIRTIGHT, and it is the user's (2026-09-08): *"if a class has two
// integrators that it needs to keep straight then SOLID::SRP dictates that class be divided."*  So the
// guarantee is not "nobody can hold two integrators" -- it is that a class holding two has a
// single-responsibility problem to fix, and the pairing is what makes that visible.
//
// ⚠ AND THE ADJOINT HALF STAYS PUBLIC ON ITS OWN.  qcMesh::MatrixOverlap has callers that want the
// adjoint ALONE and never collocate anything -- the molecular PP_Local matrix, the atom gates' 1/r and
// 1/r^2 oracles.  Making it private to this class would have broken them for a guarantee they do not
// need; a one-directional caller cannot mismatch anything.
//
// THE DENSE REALIZATION IS THE ONLY ONE HERE, deliberately.  It is the honest point-sum: no screening, no
// multigrid, no boxes -- so its two halves are adjoint by construction and it needs no proof.  The
// SCREENED realization (GPW's analytic collocation and its raw adjoint, box-truncated per multigrid level)
// implements this same interface FROM ABOVE, in the library that owns the screening data; it cannot live
// here, because qcMesh knows nothing of shells, offsets or REL_CUTOFF ladders.
module;
#include <cassert>
#include <complex>
#include <type_traits>
export module qchem.Mesh.Integrator;
export import qchem.Mesh;
export import qchem.Mesh.Quadrature;   // MatrixOverlap -- the dense adjoint
export import qchem.VectorFunction;


namespace qchem::qcMesh
{
// Module-internal scalar helpers that also work for real T -- the same pair Quadrature.C keeps
// file-private, duplicated (two lines) rather than exported: they are an implementation detail of the
// summation, not vocabulary this library wants to offer.
template <class T> inline T    IConj(const T& x) {if constexpr (std::is_floating_point_v<T>) return x; else return std::conj(x);}
template <class T> inline auto IReal(const T& x) {if constexpr (std::is_floating_point_v<T>) return x; else return std::real(x);}
}

export namespace qchem::qcMesh
{

//! \brief THE FORWARD HALF: a density matrix to its COEFFICIENTS on a fit basis.
//!
//! ★ SEGREGATED FROM THE ADJOINT (ISP, user 2026-09-09) BECAUSE THE TWO SIDES HAVE DIFFERENT CLIENTS.
//! The forward is driven by the CHARGE DENSITY -- it owns \f$D\f$ as private state and contracts it
//! against handles a caller supplies -- while the adjoint is driven by the BASIS, per block, with the
//! orbital block in hand.  Neither client wants the other's method, and a face that forces both on them is
//! the textbook ISP violation.
//!
//! ⚠ **SPLITTING THE INTERFACE DOES NOT SPLIT THE OBJECT, AND THAT IS THE WHOLE TRICK.**  The two halves
//! are two FACES of one concrete realization (\c MatrixIntegrator below), so a density holding a
//! \c MatrixForward and a basis holding a \c MatrixAdjoint that came from the SAME integrator cannot
//! mismatch -- while each names only what it uses.  Hand out two faces of one object, never two objects.
template <class T> class MatrixForward
{
public:
    virtual ~MatrixForward() = default;
    //! \brief \f$\rho_a=\sum_{ij}D_{ij}\,\langle f_a|\chi_i\chi_j\rangle=\langle f_a|\rho\rangle\f$ --
    //! the density's coefficients on my fit basis \f$\{f_a\}\f$.
    //!
    //! On a \f$\delta\f$ basis \f$f_a=\delta(r-r_a)\f$ this IS
    //! \f$\rho(r_a)=\sum_{ij}D_{ij}\overline{\chi_i(r_a)}\chi_j(r_a)\f$, which is why the two readings
    //! were never distinguished; on a plane-wave basis the same call returns \f$\tilde\rho(G_a)\f$.
    //! Real by construction for Hermitian \a D -- a density is an observable.
    virtual rvec_t Forward(const hmat_t<T>& D) const=0;

    //! \brief The SAME forward for a caller holding \f$D\f$ in FACTORED form \f$D=LL^\dagger\f$ (a thin
    //! \f$n\times r\f$ pivoted-Cholesky / natural-orbital factor):
    //! \f$\rho_g=\sum_m\left|[\Phi L]_{gm}\right|^2\f$, i.e. \f$O(n_{pts}nr)\f$ instead of
    //! \f$O(n_{pts}n^2)\f$.
    //!
    //! ★ **IT HAS A DEFAULT, AND THAT IS DELIBERATE** (user, 2026-09-09: *"any implementation of
    //! MatrixForward can support both unfactored and factored forward.  Which one (or if both) gets used
    //! should be a non-issue"*).  The default forms \f$D=LL^\dagger\f$ and calls the overload above, so
    //! **every** implementation answers both and a caller never has to ask which it has.  Overriding it is
    //! a pure OPTIMISATION -- the fast path skips materialising \f$D\f$ -- never a capability.
    //!
    //! ⚠ This replaces a runtime capability test: `applyRawFactored` is empty on realizations that cannot
    //! left-multiply a value table, and `IrrepCD::FactoredRho` asserted on it before every use.  A default
    //! makes the question disappear instead of moving it.
    virtual rvec_t Forward(const mat_t<T>& L) const
    {
        const size_t n=L.rows(), r=L.columns();
        hmat_t<T> D(n);
        for (size_t i=0;i<n;i++)
            for (size_t j=i;j<n;j++)
            {
                T s=T(0);
                for (size_t m=0;m<r;m++) s += L(i,m)*IConj(L(j,m));
                D(i,j)=s;
            }
        return Forward(D);
    }

    //! \warning ⚠ **AN IMPLEMENTATION THAT OVERRIDES ONE OVERLOAD HIDES THE OTHER.**  C++ name lookup
    //! stops at the first scope containing \c Forward, so a derived class declaring only the \c hmat_t
    //! override makes `integrator.Forward(L)` resolve to THAT one -- and since \c hmat_t is a Blaze
    //! symmetric/Hermitian ADAPTOR with a converting constructor, a thin \f$n\times r\f$ factor is then
    //! silently fed to it and Blaze throws *"Invalid setup of symmetric matrix"* at runtime, several
    //! layers from the cause.  (Measured here, 2026-09-09, the first time the pair was exercised.)
    //! ⇒ **Every implementation must say `using MatrixForward<T>::Forward;`** -- the two in this tree do,
    //! and a new one that forgets fails the `TheFactoredForwardAgreesWithTheUnfactoredOne` gate.
    virtual size_t NumCoefficients() const=0;   //!< length of the array \c Forward returns
};

//! \brief THE ADJOINT HALF: a field's coefficients back to a matrix.  See \c MatrixForward for why they
//! are separate faces of one object rather than one face or two objects.
template <class T> class MatrixAdjoint
{
public:
    virtual ~MatrixAdjoint() = default;
    //! \brief \f$\langle\chi_i|v|\chi_j\rangle=\sum_a v_a\,\langle f_a|\chi_i\chi_j\rangle\f$ for a
    //! field expanded on the SAME basis, in the same order, that the paired \c MatrixForward returns.
    //! The exact transpose of \c Forward -- same \f$\langle f_a|\chi_i\chi_j\rangle\f$, contracted the
    //! other way -- which is what makes \f$H=\partial E/\partial D\f$ hold.
    virtual hmat_t<T> Adjoint(const rvec_t& v) const=0;
    //! \f$\int f\,d^3r=\sum_a f_a\,\langle f_a|1\rangle\f$ -- the energy quadrature on the same axis.
    //! ⚠ The weights here are the FUNCTION INTEGRALS \f$\langle f_a|1\rangle\f$, not point weights; on a
    //! \f$\delta\f$ basis they coincide, which is why that distinction only surfaced when a second
    //! realization needed them (2026-09-09).
    virtual double Integrate(const rvec_t& f) const=0;
    virtual size_t NumCoefficients() const=0;   //!< length of the array \c Adjoint accepts
};

//! \brief The forward/adjoint pair of a density-matrix <-> operator assembly, as ONE object.
//!
//! Realizations differ in what they TRUNCATE (a dense point sum truncates nothing; the GPW route screens
//! per pair and boxes per multigrid level), never in semantics.  Which one a run uses is a cost decision
//! LATCHED for the run -- switching mid-SCF would change the discrete functional being minimised.
//!
//! ▶ **THIS is what a factory hands out, and its two BASES are what clients hold.**  Virtual inheritance,
//! so \c NumCoefficients is one function and the two halves provably describe the same axis.
//!
//! ★ **NOTHING HERE MENTIONS A GRID, AND THAT IS LOAD-BEARING** (user, 2026-09-09: *"if it is done right
//! (the integration grid is totally hidden inside the integrator) then this interface should also work for
//! analytic integrals"*).  What crosses these faces are COEFFICIENTS on a fit basis; whether a realization
//! reaches them by quadrature on a mesh, by an FFT, by a screened multigrid collocation, or ANALYTICALLY
//! is its own business.  That is \c doc/Pins.md pin 2 -- a fit is (integration grid) x (fit basis),
//! orthogonal axes -- with this interface parameterised on the basis alone.  \c NumPoints was renamed
//! \c NumCoefficients on 2026-09-09 for exactly this reason: it was the last word in the face that
//! presumed a grid.
template <class T> class MatrixIntegrator
    : public virtual MatrixForward<T>
    , public virtual MatrixAdjoint<T>
{
public:
    virtual ~MatrixIntegrator() = default;
};

//! \brief The DENSE realization: an honest point sum over a mesh, with no screening anywhere.
//!
//! Its two directions are adjoint BY CONSTRUCTION -- same points, same weights, same basis evaluation, no
//! truncation to get out of step -- so it is also the natural REFERENCE against which a screened
//! realization is gated.
//!
//! \warning It is \f$O(n_{pts}n^2)\f$ per direction with no sparsity, which is why production periodic runs
//! use the screened realization instead.  Do not reach for this one on a large cell because it is simple.
template <class T> class DenseMatrixIntegrator
    : public virtual MatrixIntegrator<T>
{
public:
    //! \a mesh and \a basis must outlive this object: it holds them by reference, because an integrator is
    //! a VIEW of a (mesh, basis) pair and copying either would be a lie about who owns them.
    DenseMatrixIntegrator(const Mesh& mesh, const VectorFunction<T>& basis)
        : itsMesh(mesh), itsBasis(basis) {}

    //! Un-hide the base's FACTORED overload (see the warning on \c MatrixForward::NumCoefficients): overriding
    //! the \c hmat_t one would otherwise hide it, and `Forward(L)` would convert into the wrong overload.
    using MatrixForward<T>::Forward;

    virtual rvec_t Forward(const hmat_t<T>& D) const override;
    virtual hmat_t<T> Adjoint(const rvec_t& v) const override
    {
        assert(v.size()==itsMesh.size() && "MatrixIntegrator::Adjoint: one field value per mesh point");
        return MatrixOverlap(itsMesh, itsBasis, v);
    }
    virtual double Integrate(const rvec_t& f) const override {return qcMesh::Integrate(itsMesh, f);}
    virtual size_t NumCoefficients() const override {return itsMesh.size();}

private:
    const Mesh&              itsMesh;
    const VectorFunction<T>& itsBasis;
};

//! \f$\rho(r_g)=\sum_{ij}D_{ij}\overline{\chi_i}\chi_j\f$ -- the exact adjoint of \c MatrixOverlap: the
//! same \f$\overline{\chi_i}\ldots\chi_j\f$ ordering, and the weights deliberately NOT applied here (they
//! belong to the integration, not to the density, and applying them on both sides would square them).
template <class T> rvec_t DenseMatrixIntegrator<T>::Forward(const hmat_t<T>& D) const
{
    const size_t n=itsBasis.GetVectorSize();
    assert(D.rows()==n && "MatrixIntegrator::Forward: the density matrix must match the basis");
    const rvec3vec_t& R=itsMesh.Points();
    rvec_t rho(itsMesh.size(), 0.0);
    for (size_t g=0; g<itsMesh.size(); g++)
    {
        vec_t<T> p=itsBasis(R[g]);
        T s=T(0);
        for (size_t i=0; i<n; i++)
            for (size_t j=0; j<n; j++)
                s += IConj(p[i])*D(i,j)*p[j];
        rho[g] = IReal(s);          // Hermitian D => real; the imaginary part is roundoff and is dropped
    }
    return rho;
}

} // namespace
