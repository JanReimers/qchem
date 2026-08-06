// File: BasisSet/ImplicitAngular_IBS.C  Capability: an irrep basis whose ANGULAR factor is implicit.
module;
export module qchem.BasisSet.ImplicitAngular_IBS;
import qchem.Types;   // rvec_t

export namespace qchem::BasisSet
{

//! \brief CAPABILITY: an irrep basis whose \b angular factor is IMPLICIT -- each stored function is a
//! purely RADIAL \f$\chi_i(r)\f$, and the \f$Y_{lm}\f$ that completes it to a 3-D orbital is carried by the
//! IRREP, not by the stored function.  This is the atomic lineage: a spherical solver keeps the
//! \f$m\f$-degeneracy in the OCCUPATIONS, so one radial block serves every \f$m\f$ of its \f$l\f$.
//! Consequently the block's 3-D face returns \f$\chi_i(|r|)\f$ with NO angular dependence.
//!
//! A consumer that must form \f$\langle\chi|f(r)Y_{l'm'}\rangle\f$ therefore CANNOT go through the 3-D face
//! and a 3-D quadrature: with no angular structure to integrate against, \f$\int Y_{l'm'}d\Omega\f$ makes
//! every such projection collapse to \f$\sqrt{4\pi}\int\chi f r^2dr\f$ for \f$l'=0\f$ -- spuriously coupling
//! blocks of EVERY \f$l\f$ -- and to ZERO for \f$l'\ge1\f$, silently deleting the whole projector.  It asks
//! HERE instead: match \c ImplicitL() against \f$l'\f$ (the projection vanishes unless they agree) and build
//! the radial integral from \c RadialValues.
//!
//! MEASURED consequence of getting this wrong: the occupied-d Kleinman-Bylander defect found via MnO
//! (doc/SymmetryUpgradePlan.md §7 step 7; gate \c A_PP.PerLKleinmanBylanderOracle) -- an s projector leaking
//! into the p/d/f blocks while every p/d/f projector evaluated to \f$10^{-33}\f$.
//!
//! Realized by the ATOMIC irrep basis sets only; a molecular/plane-wave block carries its angular factor
//! EXPLICITLY and simply does not implement this (the consumer keeps its 3-D path).  Reached by
//! abstract->abstract \c dynamic_cast, like the other capability faces.
class ImplicitAngular_IBS
{
public:
    virtual ~ImplicitAngular_IBS() {}

    //! The implicit angular momentum \f$l\f$ of this block (the \f$Y_{lm}\f$ factor its functions omit).
    virtual int    ImplicitL   () const=0;

    //! The block's radial functions at radius \a r: \f$\chi_i(r)\f$, one entry per basis function, in the
    //! block's own normalisation (this project: \f$\int|\chi_i|^2 d^3r=1\f$, i.e. the \f$4\pi\f$ of the full
    //! solid angle is folded IN -- so the physical radial orbital is \f$f_i=\sqrt{4\pi}\,\chi_i\f$ with
    //! \f$\int f_i f_j r^2dr=\delta_{ij}\f$).
    virtual rvec_t RadialValues(double r) const=0;
};

} //namespace
