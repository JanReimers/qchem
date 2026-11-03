// File: Symmetry/Representation.C  The matrix representation of a symmetry operation on a set of functions, and the two
// similarity transforms it induces (V-BS-SYMMETRY S3, user ruling 2026-11-03: "group-theory similarity transforms belong in
// qcSymmetry; the IBS only does them because it holds M").
//
// A set of functions {φ_I} carries a representation g -> M_g of the ops that map its span onto itself (its STABILISER):
// \f$\phi_J(g^{-1}r)=\sum_I M_{IJ}\,\phi_I(r)\f$.  Two kinds of object live on such a set and transform OPPOSITELY:
//  - MATRIX ELEMENTS \f$O_{JK}=\langle\phi_J|\hat O|\phi_K\rangle\f$ (Fock, overlap, kinetic, a potential's matrix):
//    \c Image(g,O) = \f$M_g^\dagger O M_g\f$, the matrix elements of \f$g\hat Og^{-1}\f$.
//  - COEFFICIENT-SPACE objects \f$X\f$ with \f$\rho(r)=\sum X_{IJ}\phi_I(r)\phi_J(r)^*\f$ (density matrices, transition
//    densities): \c AdjointImage(g,X) = \f$M_g X M_g^\dagger\f$, the coefficients of \f$\rho(g^{-1}r)\f$.
// They are ADJOINT under the trace pairing, \f$\mathrm{Tr}(\mathrm{AdjointImage}(g,X)\,O)=\mathrm{Tr}(X\,\mathrm{Image}(g,O))\f$,
// which is what makes an energy \f$E=\mathrm{Tr}(DF)\f$ invariant.  M need not be unitary: on a non-orthogonal set it is
// unitary in the overlap metric, \f$M^\dagger SM=S\f$.
//
// Three realizations, chosen by what M looks like (the cost argument is why there are three, not one dense one):
//  - Dense: M stored as a matrix (Gaussian AO sets: a shell permutation times small angular blocks, n_AO in the hundreds).
//  - Monomial: M a permutation times phases (plane waves: each wave goes to ONE wave), O(n^2) per op where dense is O(n^3).
//  - Identity: the set's bilinear images are the objects themselves -- a 1-D irrep block (M = χ·1, |χ|^2 = 1), an atom's l
//    block (D^l acts on the implicit m partners, 1 on the stored radial functions), or any set asked only about E.
// Every realization REFUSES an op outside the stabiliser it was built for (a throw: a broken invariant no caller can act on).
// The spin action σ of a Shubnikov op is ignored here -- it acts on the spin channel, which the caller chooses.
module;
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <cmath>
#include <complex>
export module qchem.Symmetry.Representation;
export import qchem.Symmetry.SymOp;
export import qchem.Types;
import qchem.Blaze;   // the matrix products inside the templates (ADL at instantiation needs the definition context)

export namespace qchem::Symmetry
{

//! \brief Anything that holds a representation \f$M_g\f$ of the ops in its stabiliser.  Data-free interface.
//!
//! \f$\mathrm{Tr}(\mathrm{AdjointImage}(g,X)\,O)=\mathrm{Tr}(X\,\mathrm{Image}(g,O))\f$ for every \a g in the stabiliser.
//! Both THROW (std::logic_error) for an op the set was not built for: a silent answer there would be a wrong one.
template <class T> class Representation
{
public:
    virtual ~Representation() = default;
    //! MATRIX ELEMENTS \a O (Fock, S, H): \f$M_g^\dagger O M_g\f$.
    virtual mat_t<T> Image       (const SymOp& g, const hmat_t<T>& O) const=0;
    //! COEFFICIENT-SPACE \a X (a density matrix): \f$M_g X M_g^\dagger\f$.  The result is Hermitian to roundoff; it is not
    //! symmetrized here, so a caller that sums images (a group average) rounds once, at the end.
    virtual mat_t<T> AdjointImage(const SymOp& g, const hmat_t<T>& X) const=0;
};

//! The spatial part of two ops agrees EXACTLY (σ ignored).  Op sets come from one detection, so a value key is safe.
inline bool SameSpatial(const SymOp& a, const SymOp& b)
{
    for (int i=1;i<=3;i++) for (int j=1;j<=3;j++) if (a.W(i,j)!=b.W(i,j)) return false;
    return a.tau==b.tau;
}

//! \brief WHICH ops a set's representation is defined for.  {E} by default; a shared op list (one group, many blocks); or
//! every orthogonal op fixing a point (an atom).  Shared, never copied per block.
class Stabiliser
{
public:
    Stabiliser() = default;   //!< {E}
    static Stabiliser OpList(std::shared_ptr<const std::vector<SymOp>> ops) {Stabiliser s; s.itsOps=std::move(ops); return s;}
    static Stabiliser FixingPoint(const rvec3_t& centre) {Stabiliser s; s.itsFixesPoint=true; s.itsCentre=centre; return s;}
    bool Contains(const SymOp& g) const
    {
        constexpr double tol=1e-10;
        if (itsOps) {for (const auto& op : *itsOps) if (SameSpatial(op,g)) return true; return false;}
        if (itsFixesPoint)
        {
            const Matrix3D<double> G=Transpose(g.W)*g.W;   // orthogonal (the Cartesian frame, A = 1)
            for (int i=1;i<=3;i++) for (int j=1;j<=3;j++) if (std::fabs(G(i,j)-(i==j?1.0:0.0))>tol) return false;
            const rvec3_t d=g.W*itsCentre+g.tau-itsCentre;
            return std::fabs(d.x)+std::fabs(d.y)+std::fabs(d.z)<=tol;
        }
        for (int i=1;i<=3;i++) for (int j=1;j<=3;j++) if (std::fabs(g.W(i,j)-(i==j?1.0:0.0))>tol) return false;
        return std::fabs(g.tau.x)+std::fabs(g.tau.y)+std::fabs(g.tau.z)<=tol;
    }
private:
    std::shared_ptr<const std::vector<SymOp>> itsOps;
    bool    itsFixesPoint=false;
    rvec3_t itsCentre=rvec3_t(0,0,0);
};

//! \brief The set's bilinear images ARE the objects: Image(g,O)=O, AdjointImage(g,X)=X, for \a g in the stabiliser.
//! Exact for a 1-D irrep block (\f$M_g=\chi(g)\,1\f$ and \f$|\chi|^2=1\f$, real or complex characters alike), for an atom's l
//! block (\f$D^l(g)\f$ acts on the implicit m partners, 1 on the stored radial functions; the \f$(-1)^l\f$ and \f$D^l\f$
//! phases cancel in radial-space bilinears), and for any set asked only about E.  Stores no M; THROWS for any other op.
template <class T> class IdentityRepresentation : public virtual Representation<T>
{
public:
    explicit IdentityRepresentation(Stabiliser s={}) : itsStabiliser(std::move(s)) {}
    virtual mat_t<T> Image       (const SymOp& g, const hmat_t<T>& O) const override {Require(g); return mat_t<T>(O);}
    virtual mat_t<T> AdjointImage(const SymOp& g, const hmat_t<T>& X) const override {Require(g); return mat_t<T>(X);}
private:
    void Require(const SymOp& g) const
    {
        if (!itsStabiliser.Contains(g))
            throw std::logic_error("IdentityRepresentation: asked about an op outside the stabiliser this set was built for -- "
                                   "its functions move under it, and answering M = 1 would be wrong");
    }
    Stabiliser itsStabiliser;
};

//! \brief M stored as a matrix, built by the set on first use per op (\c MakeRepMatrix, the pull hook) and cached by the op's
//! value.  Building M is the set's business (it knows its functions); caching and applying it are this class's.
template <class T> class DenseRepresentation : public virtual Representation<T>
{
public:
    virtual mat_t<T> Image       (const SymOp& g, const hmat_t<T>& O) const override
    {
        const mat_t<T>& M=RepMatrix(g);
        return blazem::ctrans(M)*mat_t<T>(O)*M;
    }
    virtual mat_t<T> AdjointImage(const SymOp& g, const hmat_t<T>& X) const override
    {
        const mat_t<T>& M=RepMatrix(g);
        return M*mat_t<T>(X)*blazem::ctrans(M);
    }
protected:
    //! \f$M_g\f$ on this set's functions; THROWS if \a g does not map the set onto itself.
    virtual mat_t<T> MakeRepMatrix(const SymOp& g) const=0;
private:
    const mat_t<T>& RepMatrix(const SymOp& g) const
    {
        for (const auto& [op,M] : itsCache) if (SameSpatial(op,g)) return M;
        itsCache.emplace_back(g, MakeRepMatrix(g));
        return itsCache.back().second;
    }
    mutable std::vector<std::pair<SymOp, mat_t<T>>> itsCache;   //!< per op, by value; empty until an op is asked about
};

//! \brief M a permutation times phases, \f$(M_g)_{\pi(i),i}=p_i\f$ -- a plane-wave set, where each wave goes to ONE wave.  The
//! images are re-indexings, O(n²) per op; the dense product would be O(n³) and n² of storage per op.  Complex only.
class MonomialRepresentation : public virtual Representation<dcmplx>
{
public:
    //! One op's rep: function i goes to function perm[i] with amplitude phase[i].
    struct Monomial { std::vector<size_t> perm; cvec_t phase; };
    virtual mat_t<dcmplx> Image       (const SymOp& g, const hmat_t<dcmplx>& O) const override
    {
        const Monomial& m=GetMonomial(g);
        const size_t n=m.perm.size();
        mat_t<dcmplx> R(n,n);
        for (size_t i=0;i<n;i++) for (size_t j=0;j<n;j++) R(i,j)=std::conj(m.phase[i])*m.phase[j]*O(m.perm[i],m.perm[j]);
        return R;
    }
    virtual mat_t<dcmplx> AdjointImage(const SymOp& g, const hmat_t<dcmplx>& X) const override
    {
        const Monomial& m=GetMonomial(g);
        const size_t n=m.perm.size();
        mat_t<dcmplx> R(n,n);
        for (size_t i=0;i<n;i++)
            for (size_t j=0;j<n;j++)
            {
                const dcmplx f=m.phase[i]*std::conj(m.phase[j]);
                R(m.perm[i],m.perm[j])=f*X(i,j);
            }
        return R;
    }
protected:
    //! \f$M_g\f$ as a monomial on this set's functions; THROWS if \a g does not map the set onto itself.
    virtual Monomial MakeMonomial(const SymOp& g) const=0;
private:
    const Monomial& GetMonomial(const SymOp& g) const
    {
        for (const auto& [op,m] : itsCache) if (SameSpatial(op,g)) return m;
        itsCache.emplace_back(g, MakeMonomial(g));
        return itsCache.back().second;
    }
    mutable std::vector<std::pair<SymOp, Monomial>> itsCache;
};

} // namespace qchem::Symmetry
