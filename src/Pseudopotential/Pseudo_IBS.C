// File: Pseudopotential/Pseudo_IBS.C  Abstract reciprocal-space PSEUDOPOTENTIAL assembly capability.
//
// The external one-body potential of a pseudopotential crystal, V_ext = V_loc + V_NL, assembled in
// reciprocal (G-)space.  This is the basis-side half of the pseudo-wall: the TERM (Hamiltonian-side)
// owns the pseudopotential MODEL and hands it here; the basis owns the G-space assembly -- the
// structure factor e^{-iG.tau}, the 1/Omega, the G=0 handling.  Lifted OFF Band_FT_IBS so qcBasisSet
// names no pseudopotential type at all; this capability lives with the models, in qcPseudopotential,
// and a plane-wave basis (qcLattice_BS) realizes it.  A term reaches it the sanctioned way: holding the
// abstract orbital basis and dynamic_cast-ing ACROSS to this capability (both abstract).
module;
export module qchem.Pseudopotential.Pseudo_IBS;
export import qchem.Pseudopotential.LocalPotential;
export import qchem.Pseudopotential.SeparablePotential;
import qchem.Structure;   // Structure -- the external-potential structure-factor source (positions, Z)
import qchem.Types;       // hmat_t<dcmplx>

export namespace Pseudopotential
{

//! \brief A basis that assembles the EXTERNAL (pseudo)potential matrices in reciprocal space from a
//! pseudopotential model.  Only the complex (plane-wave) path realizes this; the basis owns the
//! G-space assembly, the held model supplies the per-species physics.  \f$V_{ext}=V_{loc}+V_{NL}\f$.
class Pseudo_IBS
{
public:
    //! \brief Local external-potential matrix \f$\langle i|V_{loc}|j\rangle\f$ for \a structure, from the
    //! local model \a loc: \f$\frac1\Omega\sum_a v_{loc}(Z_a,|\Delta G|^2)\,e^{-i\Delta G\cdot\tau_a}\f$,
    //! \f$\Delta G\ne0\f$ (the \f$\Delta G=0\f$ term is the dropped neutralising-background shift).  Hermitian.
    virtual hmat_t<dcmplx> MakeLocalPotential(const Structure*, const LocalPotential& loc) const=0;

    //! \brief Separable (Kleinman-Bylander) nonlocal matrix \f$\langle i|V_{NL}|j\rangle\f$ for \a structure
    //! from the projector model \a nl: \f$\frac1\Omega\sum_a e^{-i\Delta G\cdot\tau_a}\sum_p\tilde\beta_p(|k+G|)
    //! D_p\tilde\beta_p(|k+G'|)\f$ (with the \f$(2l+1)P_l\f$ angular weight per channel).  Hermitian.
    virtual hmat_t<dcmplx> MakeSeparablePotential(const Structure*, const SeparablePotential& nl) const=0;

    //! \brief Energy carried by the local potential's DROPPED \f$G=0\f$ component for the model \a loc and a
    //! density of \a numElectrons electrons: the uniform electron-ion alignment \f$(N/\Omega)\sum_a\alpha_a\f$,
    //! \f$\alpha_a=\int[V_{loc}^a+Z_a/r]\f$ (the basis supplies \f$\Omega\f$; the model supplies \f$\alpha\f$).
    //! Enters the total energy but NOT the Hamiltonian matrix (the \f$G=0\f$ potential is a dropped constant).
    virtual double ExternalG0Energy(const Structure*, const LocalPotential& loc, double numElectrons) const=0;
};

}//namespace
