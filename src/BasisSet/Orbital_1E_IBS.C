// File: BasisSet/Orbital_1E_IBS.C Orbital that knows enough integrals for a 1 electron calculation.
module;
export module qchem.BasisSet.Orbital_1E_IBS;
export import qchem.BasisSet.IrrepBasisSet;
export import qchem.Structure;

export namespace qchem::BasisSet
{


//--------------------------------------------------------------------------------
//
//! The method of integral evaluation is of course strongly dependant on the
//! precise details of basis functions or basis set.  
//! All integral functions except MakeNormalization return normalized integrals.
//! Interfaces for 1 electron integrals used for all Irrep basis sets: 1E,Fit,HF,DFT,DHF  
//! The calls return matrix refrences which implies they are buffered behind the scenes.
//

//! \brief Interface for the kinetic-energy BUILDING BLOCK -- which is NOT the kinetic energy.
//!
//! MEANING OF "Kinetic" here:  MakeKinetic()/Kinetic() return the matrix of \f$-\nabla^2\f$, i.e.
//! \f[ \langle p^2\rangle \;=\; \langle a|-\nabla^2|b\rangle \;=\; \int d^3\vec r\,\nabla g_a\cdot\nabla g_b \f]
//! the positive-definite "Grad2" / momentum-squared block.  It does \b NOT include the \f$\tfrac12\f$.
//! (The previous doc here claimed \f$\langle a|-\tfrac12\nabla^2|b\rangle\f$ -- that was wrong and is
//! the likely origin of much confusion: the \f$\tfrac12\f$ lives at the energy boundary, not here.)
//!
//!  - Non-relativistic kinetic ENERGY \f$T=-\tfrac12\nabla^2\f$: the \f$\tfrac12\f$ is applied in
//!    Hamiltonian/Internal/Imp/Kinetic.C as \f$0.5\cdot\f$Kinetic().
//!  - Dirac: DiracKinetic.C / the RKB chain consume this \f$\langle p^2\rangle\f$ block with their
//!    own \f$c_{light}\f$ factors (no \f$\tfrac12\f$) -- those factors are unverified (see DHF below).
//!
//! (A later cleanup will rename this block to P2()/MakeP2() so the name stops saying "Kinetic".)
template <class T> class Integrals_Kinetic : public virtual IrrepBasisSet_IDs
{
public:
    //! \return \f$\langle p^2\rangle=\langle -\nabla^2\rangle\f$ (NO 1/2). Built once per {radial,angular} ID pair.
    virtual hmat_t<T>  MakeKinetic() const=0;
    //! Cached accessor for the \f$\langle p^2\rangle\f$ block (see MakeKinetic()).
    const   hmat_t<T>&     Kinetic() const;

    //! \name d/dR (the user's ruling 2026-11-08: a derivative version beside a hub face's method; GeometryOptimizationPlan §13.2)
    //!@{
    //! Can this block give MakeKineticPullback?  CLIENT: the kinetic term's MissingGradients() entry at setup.  False exactly when it would refuse.
    virtual bool CanPullbackKinetic() const=0;
    //! THE PULLBACK of K = <-∇²> (no ½, as MakeKinetic): Σ_ij Re(D_ij* ∂K_ij/∂R_A) for every atom A of \a st, in its site order.  Exact zeros
    //! for functions that do not move with the atoms (plane waves); NoGradient for a family with no derivative yet.
    virtual Outcome<rvec3vec_t,AD::NoGradient> MakeKineticPullback(const hmat_t<T>& D, const Structure* st) const=0;
    //! THE STRAIN DERIVATIVE beside it (as Integrals_Overlap::MakeOverlapVirial; CellHalf_Design §5, user ruling §6.2): V_ab = ∂Re Tr(D† K)/∂ε_ab
    //! at fixed D and fractional k (no ½), the channel-A image virial for functions that ride the atoms.  CLIENT: the kinetic term's virial.
    //! NOT pure: a family without a strain derivative refuses by default.
    virtual bool CanVirialKinetic() const {return false;}
    virtual Outcome<rmat3d_t,AD::NoGradient> MakeKineticVirial(const hmat_t<T>&, const Structure*) const
        {return NoGeometryDerivativeYet<rmat3d_t>("this basis block", "kinetic virial (strain derivative)");}
    //!@}
};

//! \brief Interface for electron-nucleus attraction integrals.
//! Nuclear attraction \f$ \sum_{i}\left\langle a\left|\frac{-Z_{i}}{\left|\vec{r}-\vec{R}_{c}\right|}\right|b\right\rangle =-\sum_{i}Z_{i}\int d^{3}\vec{r}\:g_{a}\left(\vec{r}\right)\frac{1}{\left|\vec{r}-\vec{R}_{c}\right|}g_{b}\left(\vec{r}\right)\f$
template <class T> class Integrals_Nuclear : public virtual IrrepBasisSet_IDs
{
public:
    virtual hmat_t<T>  MakeNuclear(const Structure* cl) const=0; //Only called once for a given {radial,angular, structure} ID triple.
    const   hmat_t<T>&     Nuclear(const Structure* cl) const;

    //! \name d/dR (the user's ruling 2026-11-08; GeometryOptimizationPlan §13.2)
    //!@{
    //! Can this block give MakeNuclearPullback?  CLIENT: the nuclear-attraction term's MissingGradients() entry at setup.
    virtual bool CanPullbackNuclear() const=0;
    //! THE PULLBACK of V: Σ_ij Re(D_ij* ∂V_ij/∂R_A) for every atom A of \a cl, with the basis centres AND the nuclei moving.  On plane waves
    //! the basis does not move but V does: the structure-factor Hellmann-Feynman force.  NoGradient for a family with no derivative yet.
    virtual Outcome<rvec3vec_t,AD::NoGradient> MakeNuclearPullback(const hmat_t<T>& D, const Structure* cl) const=0;
    //!@}
};


//
// Define an orbital irrep basis set which supports integrals for 1-electron (1E) orbital calculations.
// Mix-in the integral interfaces required for a 1E orbital basis. 
//
template <class T> class Orbital_1E_IBS
    : public virtual Evaluatable_IBS<T> //brings in symmetry + op()(r): an orbital basis IS evaluatable
    , public virtual Integrals_Overlap<T>
    , public virtual Integrals_Kinetic<T>
    , public virtual Integrals_Nuclear<T>
{
public:
};

typedef Orbital_1E_IBS<double>    Real_OIBS;
typedef Orbital_1E_IBS<dcmplx> Complex_OIBS;

} //namespace
