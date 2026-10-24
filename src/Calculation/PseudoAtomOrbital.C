// File: Calculation/PseudoAtomOrbital.C  THE PSEUDO-ATOM'S OWN l ORBITAL -- what a +U projector is
// (doc/HubbardProjectorPlan.md D5).
//
// QE and hp.x define a Hubbard projector as an atomic function (the UPF's PP_CHI, which our own gth2upf wrote
// from this same GTH pseudo-atom).  So the projector is built here, once per species, in the pseudo-atom's OWN
// exponents -- a complete 16-exponent even-tempered pool, LDA, unpolarized -- and handed to the +U term as a
// closed Gaussian radial.  Nothing here knows the solid's orbital basis: that independence is the point (the
// site-shell projection it replaces changed U with the basis, doc/Records/HubbardUHistory.md 2026-10-24).
//
// The pool reproduces CP2K's ATOM code to 3 / 0.3 / 0.1 mHa on Mn (E, eps 4s, eps 3d;
// UTCalculation ValenceBasisGen.MnQ7PseudoAtomInALargePool).
module;
#include <string>
#include <vector>
export module qchem.PseudoAtomOrbital;
export import qchem.Math.Gaussian;   // Math::Gaussian (the radial a +U manifold carries)

export namespace qchem
{
//! \brief The lowest occupied \f$l\f$ orbital of the GTH pseudo-atom (LDA, unpolarized) in a complete
//! even-tempered pool: \f$R(r)=\sum_kc_kN_kr^le^{-\alpha_kr^2}\f$ over unit-normalised primitives (\f$N_k\f$).
struct PseudoAtomOrbital
{
    std::string         element;
    int                 Zion = 0;
    int                 l    = 0;
    std::vector<double> exponents;      //!< the pool \f$\alpha_k\f$
    std::vector<double> coefficients;   //!< \f$c_k\f$ over unit-normalised primitives (\f$\sum c_kc_k'S_{kk'}=1\f$)
    double              eps  = 0.0;     //!< the orbital's eigenvalue (Ha)
    double              Eatom= 0.0;     //!< the pseudo-atom's total energy (Ha)
    //! The radial as raw Gaussian terms \f$\{c_kN_k,\,n=0,\,\alpha_k\}\f$, \f$r^l\f$ carried by the owner -- the
    //! form \c HubbardManifold::atomicFunction and \c BasisSet::AtomCentredSet_Gaussian take.
    std::vector<Math::Gaussian> Radial() const;
};

//! \brief The pseudo-atom orbital of \a element at valence \a Zion (0 = the GTH default) and angular momentum
//! \a l, cached per (element, Zion, l, pool).  \a lmaxPool: the highest l the pool carries; < 0 = the default,
//! \f$\max(l,2)\f$ (\f$\max(l,3)\f$ from La on) -- every channel an occupied valence shell can need; an empty
//! channel changes nothing in a spherical unpolarized atom.  THROWS if the atom does not converge or has no
//! occupied \f$l\f$ orbital.
const PseudoAtomOrbital& PseudoAtomOrbitalOf(const std::string& element, int Zion, int l, int lmaxPool=-1);

} //namespace
