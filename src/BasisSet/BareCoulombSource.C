// File: BasisSet/BareCoulombSource.C  "I can hand you the BARE two-electron integrals over a chosen subset of
// my functions": the capability face the ACBN0 estimator of DFT+U needs (doc/OpenWork.md step 5 increment 3,
// 2026-09-21) -- the on-site Hartree-Fock energy of a Hubbard manifold is built from (m1 m2|m3 m4) over the
// manifold's own functions, unscreened and in the CENTRAL CELL ONLY (Agapito et al. 2015 eq 4 with g=l=m=0):
// the renormalised density matrix supplies the screening, so the integrals must NOT be lattice-summed.
//
// WHY A FACE OF ITS OWN.  The four-centre machinery the Fock build uses (Orbital_HF_IBS) delivers whole
// symmetry-packed ERI4 blocks for a Coulomb/exchange CONTRACTION; a client wanting the raw m^4 tensor over a
// handful of functions is a different integral CONSUMER, and the pseudo-wall rule (change a basis interface
// only for a NEW INTEGRAL TYPE) is exactly what this is.  Answered where the integrals are analytic (the
// Gaussian FourC mixin), FORWARDED by a Bloch block to its molecular block (a Bloch sum of an AO is the AO in
// the central cell), TRANSFORMED by a spherical view through its own cart->sphere map.  A block that cannot
// answer does not derive the face; a consumer cross-casts abstract->abstract and throws with a clear message.
module;
#include <cstddef>
#include <vector>
export module qchem.BasisSet.BareCoulombSource;
import qchem.Types;

export namespace qchem::BasisSet
{

//! \brief The \f$m^4\f$ tensor \f$(m_1m_2|m_3m_4)=\int\phi_{m_1}\phi_{m_2}\,r_{12}^{-1}\,\phi_{m_3}\phi_{m_4}\f$ over
//! \f$m\f$ functions, dense and row-major -- small by construction (a manifold), so the 8-fold symmetry is not
//! packed away; the client's contractions read it in every index order.
class ERI4Block
{
public:
    ERI4Block() : itsM(0) {}
    explicit ERI4Block(size_t m) : itsM(m), itsV(m*m*m*m, 0.0) {}
    size_t Size() const {return itsM;}
    double  operator()(size_t a, size_t b, size_t c, size_t d) const {return itsV[((a*itsM+b)*itsM+c)*itsM+d];}
    double& operator()(size_t a, size_t b, size_t c, size_t d)       {return itsV[((a*itsM+b)*itsM+c)*itsM+d];}
    const rvec_t& Data() const {return itsV;}
    //! \f$(m_1m_2|m_3m_4)\f$ carried through a linear map \f$\phi'_n=\sum_a T(a,n)\phi_a\f$ on ALL FOUR indices
    //! (a spherical view's cart->sphere map, a contraction, a rotation): \f$T\f$ is \f$m\times m'\f$.
    ERI4Block Transform(const rmat_t& T) const;
private:
    size_t itsM;
    rvec_t itsV;
};

class BareCoulombSource
{
public:
    virtual ~BareCoulombSource() = default;
    //! \brief The bare integrals over the functions \a cols of THIS block (its own function indices), in
    //! \a cols order.  No lattice sum, no screening, no metric: the functions as they are in the central cell.
    virtual ERI4Block BareCoulomb(const std::vector<size_t>& cols) const = 0;
};

} // namespace
