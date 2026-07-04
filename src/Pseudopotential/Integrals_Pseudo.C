// File: Pseudopotential/Integrals_Pseudo.C  Abstract PSEUDOPOTENTIAL operator-assembly capability.
//
// The external one-body potential of a pseudopotential, V_ext = V_loc + V_NL, assembled into its matrix
// blocks in a basis.  A capability MIXIN (cf. Integrals_Overlap/Kinetic/Nuclear) -- NOT itself a basis
// set: a basis inherits it to gain pseudopotential support.  This is the basis-side half of the
// pseudo-wall: the TERM (Hamiltonian-side) owns the pseudopotential MODEL and hands it here; the basis
// owns the assembly (for plane waves: the structure factor e^{-iG.tau}, the 1/Omega, the G=0 handling).
// It lives with the models in qcPseudopotential so qcBasisSet names no pseudopotential type; the plane-
// wave basis realizes Integrals_Pseudo<dcmplx>, and a future molecular ECP would realize the <double>
// instantiation.  A term reaches it the sanctioned way: holding the abstract orbital basis and
// dynamic_cast-ing ACROSS to this capability (both abstract).
module;
export module qchem.Pseudopotential.Integrals_Pseudo;
export import qchem.Pseudopotential.LocalPotential;
export import qchem.Pseudopotential.SeparablePotential;
import qchem.Structure;   // Structure -- the external-potential structure-factor source (positions, Z)
import qchem.Types;       // hmat_t<T>

export namespace qchem::Pseudopotential
{

//! \brief Pseudopotential operator-assembly mixin: builds the external (pseudo)potential matrix blocks
//! from a pseudopotential model.  \tparam T the matrix scalar (dcmplx for a plane-wave basis; double for
//! a future real-Gaussian molecular ECP).  The basis owns the assembly; the held model (term-side)
//! supplies the per-species physics.  \f$V_{ext}=V_{loc}+V_{NL}\f$ (the separable part is optional -- a
//! model fact, nl absent -- not a basis-capability split, so local + separable live in one interface).
template <class T> class Integrals_Pseudo
{
public:
    //! \brief Local external-potential matrix \f$\langle i|V_{loc}|j\rangle\f$ for \a structure, from the
    //! local model \a loc: \f$\frac1\Omega\sum_a v_{loc}(Z_a,|\Delta G|^2)\,e^{-i\Delta G\cdot\tau_a}\f$,
    //! \f$\Delta G\ne0\f$ (the \f$\Delta G=0\f$ term is the dropped neutralising-background shift).  Hermitian.
    virtual hmat_t<T> MakeLocalPotential(const Structure*, const LocalPotential& loc) const=0;

    //! \brief Separable (Kleinman-Bylander) nonlocal matrix \f$\langle i|V_{NL}|j\rangle\f$ for \a structure
    //! from the projector model \a nl: \f$\frac1\Omega\sum_a e^{-i\Delta G\cdot\tau_a}\sum_p\tilde\beta_p(|k+G|)
    //! D_p\tilde\beta_p(|k+G'|)\f$ (with the \f$(2l+1)P_l\f$ angular weight per channel).  Hermitian.
    virtual hmat_t<T> MakeSeparablePotential(const Structure*, const SeparablePotential& nl) const=0;
};

}//namespace
