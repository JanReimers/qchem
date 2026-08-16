// File: Orbital.C  Interface for Orbitals.
module;
#include <vector>
#include <memory>
export module qchem.Orbitals;
export import qchem.ChargeDensity;
export import qchem.Symmetry.Orbital;

import qchem.ScalarFunction;
import qchem.VectorFunction;
import qchem.Streamable;
import qchem.Iterators;

export namespace qchem::Orbitals
{

using qchem::ChargeDensity::rDM_CD;
using qchem::ChargeDensity::tDM_CD;

//#############################################################
//
//  An orbital is 1) a real space function, 2) it holds some number
//  and type of electrons, 3) it has an eigen energy.  This non
//  templated portion is independant of whether the orbital real
//  or complex valued.
//
class Orbital
    : public virtual Streamable
{
public:
    virtual ~Orbital() {};
    virtual bool   IsOccupied   (       ) const=0;
    virtual double GetOccupation(       ) const=0;
    virtual void   Empty        (       )      =0;
    virtual double TakeElectrons(double )      =0;
    //! Set this orbital's occupation directly (the whole level, i.e. degeneracy*f for a Fermi-smeared
    //! state).  Aufbau uses TakeElectrons; Fermi smearing sets g*f_i explicitly (doc/GPWPlan1.md 4b).
    virtual void   SetOccupation(double )      =0;
    virtual int    GetDegeneracy(       ) const=0;
    
    virtual double      GetEigenEnergy () const=0;
    virtual Orbital_QNs GetQNs         () const=0; //Should have principle QN + spin QN + any symmetry QNs.
    virtual std::string GetLabel       () const=0; //A text version of the QNs.
};

//---------------------------------------------------------
//
//  Templated depending or whether it is a real or
//  complex valued orbital.
//
template <class T> class TOrbital
    : public virtual Orbital
    , public virtual ScalarFunction<T>
{
public:
    virtual void AddDensityMatrix(hmat_t<T> & D, hmat_t<T> & DPrime) const=0;  // DM is Hermitian (=symmetric for real T)
    //! Coefficients in the *orthonormal* basis (C'); the metric there is the identity, so MOM
    //! orbital overlaps are plain dot products of these vectors.
    virtual const vec_t<T>& GetCoeffPrime() const=0;
};



//---------------------------------------------------------------------------
//
//  A group of orbitals is usually for one irreducable representation.
//  The most interesting member function is GetChargeDensity().  This non
//  templated portion is independant of whether the orbitals are real
//  or complex valued.
//
class Orbitals : public virtual Streamable
{
public:
    virtual ~Orbitals() {};
    virtual size_t         GetNumOrbitals     (               ) const=0;
    virtual size_t         GetNumOccOrbitals  (               ) const=0;
    virtual double         GetEigenValueChange(const Orbitals&) const=0;
    // GetChargeDensity moved to TOrbitals<T> (it returns the T-typed density tDM_CD<T>*; a complex
    // orbital group yields a complex density, which the non-templated base cannot express).
    //! This will hold spin and symmetry QNs, without the principle QN.
    virtual Irrep      GetQNs() const=0;

    //  Iterate() with no type argument yields the base Orbital* directly (no
    //  cast); Iterate<D>() dynamic_cast's each to the requested derived type D.
    //  Both come in const and mutable flavours (the latter for e.g. Empty()),
    //  over storage kept private by the concrete Orbitals.
    auto Iterate() const {return IndexProxy<const Orbitals>(this, GetNumOrbitals());}
    auto Iterate()       {return IndexProxy<      Orbitals>(this, GetNumOrbitals());}
    template <class D> auto Iterate() const {return D_IndexProxy<const D, const Orbitals>(this, GetNumOrbitals());}
    template <class D> auto Iterate()       {return D_IndexProxy<D, Orbitals>(this, GetNumOrbitals());}

    const Orbital* operator[](size_t i) const {return GetOrbital(i);}
          Orbital* operator[](size_t i)       {return GetOrbital(i);}

protected:
    virtual const Orbital* GetOrbital(size_t) const=0; //the only storage-specific
    virtual       Orbital* GetOrbital(size_t)      =0; //primitives
};

//---------------------------------------------------------
//
//  Templated depending or whether it is a real or
//  complex valued orbital.
//
//! Fermi–Dirac occupancy fraction \f$f=1/(1+e^{(e-\mu)/kT})\f$, overflow-guarded (\f$|x|>40\f$ saturates to
//! 0/1).  \a kT>0.  The ONE shared occupancy kernel for both the per-block μ-solve (\c TakeElectronsFermi) and
//! the composite cross-k global-μ fill (doc/GPWPlan1.md items 3 / 4b) -- keep them bit-identical here.
double FermiOccupancy(double e, double mu, double kT);
//! Solve the chemical potential μ by bisection so \f$\sum_i g_i\,\mathrm{FermiOccupancy}(e_i,\mu,kT)=\text{target}\f$.
//! The energies \a e are EFFECTIVE (bare ε plus any MOM eShift, already folded in by the caller); \a g are the
//! per-level capacities -- the degeneracy for a single block, or BZ-weight×degeneracy for the composite
//! cross-k fill, so ONE unweighted solver serves both (the global weight rides in \a g).  Requires
//! \a e.size()==g.size()>0 and \f$\sum g\ge\text{target}\f$ (capacity).  Returns μ (doc/GPWPlan1.md item 3).
double FermiLevel(const rvec_t& e, const rvec_t& g, double target, double kT);

template <class T> class TOrbitals
    : public virtual Orbitals
    , public virtual VectorFunction<T>
{
public:
    virtual size_t GetVectorSize() const
    {
        return GetNumOrbitals();
    }
    typedef std::tuple<double,hmat_t<T>> ds_t;   // {remaining electrons, DPrime} (DPrime Hermitian)

    virtual void  UpdateOrbitals(const mat_t<T>& U, const mat_t<T>& UPrime, const rvec_t& e)=0;
    virtual ds_t TakeElectrons (double ne)=0;
    //! MOM occupation: fill the highest-\a priority orbitals first (one score per orbital, stored order)
    //! instead of lowest-energy -- occupied-subspace continuity for a within-irrep level crossing
    //! (doc/GPWPlan §0b″).  Ties keep the stored (energy) order.  \a priority.size()==GetNumOrbitals().
    virtual ds_t TakeElectrons (double ne, const rvec_t& priority)=0;
    //! Fermi-Dirac occupation (doc/GPWPlan1.md 4b): solve the chemical potential μ by bisection so
    //! \f$\sum_i g_i f_i(\varepsilon_i)=n_e\f$ with \f$f_i=1/(1+e^{(\varepsilon_i-\mu)/kT})\f$, set each
    //! orbital's occupation to \f$g_i f_i\f$, and build D/D'.  \a kT>0 (Hartree).  Returns
    //! {\a MinusTS, DPrime} where \a MinusTS \f$=kT\sum_i g_i[f_i\ln f_i+(1-f_i)\ln(1-f_i)]\le0\f$ is the
    //! Mermin free-energy term \f$-TS\f$ (reuses the ds_t tuple; here the double is −TS, NOT leftover
    //! electrons -- μ is solved so there is no leftover).  The fractional occupations flow through the
    //! existing density build unchanged (AddDensityMatrix already scales |C⟩⟨C| by the occupation).
    virtual ds_t TakeElectronsFermi (double ne, double kT)=0;
    //! \brief The chemical potential of the LAST Fermi fill -- solved by \c TakeElectronsFermi (this block's
    //! own μ, one per constrained channel) or handed in by \c SetFermiOccupationsAtMu (a shared reservoir's).
    //! NaN before any Fermi fill has run.
    //!
    //! It is REPORTED, not just used, because Δμ between channels is the diagnostic that says whether a
    //! fixed-(nUp,nDn) fill is doing work: two μ ARE the Lagrange multipliers of the two count constraints,
    //! so their difference is the field holding the moment apart, and a run that never prints them cannot
    //! tell a free moment from a constrained one.  (Had this been printed, MnO's 27 mHa would have been read
    //! off a trace instead of inverted by hand out of an eigen table -- user, 2026-08-10.)
    //! \note Only sharp where the occupations are FRACTIONAL.  With integer fills μ is pinned no better than
    //! the channel's HOMO-LUMO gap, so compare gaps, not digits.
    virtual double GetChemicalPotential() const=0;
    //! MOM-masked Fermi (doc/GPWPlan1.md 4b, MOM+smearing composition): Fermi-fill on EFFECTIVE energies
    //! ε_i + eShift_i instead of the bare ε_i.  \a eShift (one per orbital, stored order; size 0 == all-zero
    //! == the plain overload above) carries a MOM-overlap penalty that pushes LOW-overlap states (diffuse
    //! ghosts) UP so they stay empty by CHARACTER, while the retained high-overlap states smear by their
    //! TRUE energy -- keeping the ghost out (like hard MOM) yet making occupation continuous at the frontier.
    virtual ds_t TakeElectronsFermi (double ne, double kT, const rvec_t& eShift)=0;
    //! Set Fermi occupations at a GIVEN chemical potential μ (no per-block solve): occ_i = g_i·FermiOccupancy(
    //! ε_i+eShift_i, μ, kT), accumulate this block's −TS, and build D/D'.  The composite cross-k global-μ fill
    //! solves ONE μ across the mesh (via \c FermiLevel over the aggregate) then calls this on EACH block so the
    //! charge sloshes between k-points under a single Fermi level (doc/GPWPlan1.md item 3).  \c TakeElectronsFermi
    //! is exactly {solve μ for this block's ne} then this.  \a kT>0; \a eShift size 0 == all-zero (plain Fermi).
    virtual ds_t SetFermiOccupationsAtMu (double mu, double kT, const rvec_t& eShift)=0;
    //! \brief BUILDS this orbital set's density matrix -- it ALLOCATES, hence the owning return (V1.25).
    //! It was a raw `tDM_CD<T>*` from a `Get*`, which reads as an accessor handing back a member; every
    //! caller had to know otherwise, and AtomCalculation::TotalCharge did not, leaking a whole composite
    //! per call.  Ownership now lives in the TYPE.
    virtual std::unique_ptr<tDM_CD<T>> GetChargeDensity() const=0;

    //! Per-basis-function occupation-weighted Mulliken gross population \f$P_i=(DS)_{ii}\f$, given this irrep's
    //! AO overlap \a S (supplied by the caller's IrrepBasisSet).  The density \f$D=\sum_{occ} n|C\rangle\langle
    //! C|\f$ already folds in the occupations, so \f$\sum_i P_i=\mathrm{Tr}(DS)\f$ = the irrep's electron count.
    //! A basis-USAGE diagnostic -- which primitives actually carry the density -- for tuning a valence window
    //! (peak at a window boundary => extend it; ~0 there => shrink it).  Keeps \f$D\f$ hidden: the caller hands
    //! in only its own \a S.  Real part taken (exact for a real basis; the physical gross population otherwise).
    virtual rvec_t GetBasisPopulations(const hmat_t<T>& S) const=0;

};

} //namespace