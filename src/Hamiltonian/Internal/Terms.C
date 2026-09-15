// Hamiltonian/Internal/Terms.C  Declare and export all Hamiltonian term types.
module;
#include <iosfwd>
#include <memory>
#include <map>
#include <string>
#include <vector>
#include <functional>
export module qchem.Hamiltonian.Internal.Terms;
import qchem.Hamiltonian.Internal.Term;
import qchem.Hamiltonian.Internal.ExFunctional;
import qchem.Structure;
import qchem.Fitting.FunctionFitter;          // Fitting::FunctionFitter (composed; clients never see the impl)
import qchem.ChargeDensity;
import qchem.FittedCD;
import qchem.Hamiltonian.Types;
import qchem.BasisSet.Orbital_PP_IBS;          // SpeciesRadialField / SpeciesProjectorSet_R (the real-space views the molecular PP terms quadrature)
import qchem.Mesh;                             // qcMesh::MeshParams (the quadrature mesh spec)
import qchem.BasisSet.ImplicitAngular_IBS;      // the radial/implicit-Y_lm capability (atomic KB route)


export namespace qchem::Hamiltonian
{

using ChargeDensity::ChannelOf;      // the spin channels of a density, through the face (V1.37)
using ChargeDensity::DM_ChannelOf;

// The non-relativistic kinetic ENERGY term is now the T-templated Kinetic<T>
// (qchem.Hamiltonian.Internal.Kinetic); the molecular Hamiltonians build Kinetic<double>.

// Relativistic kinetic ENERGY term (Dirac \f$c\,\vec\sigma\cdot\vec p\f$); consumes the RKB-assembled
// relativistic kinetic block directly (no 1/2). See Imp/DiracKinetic.C for the unverified-factor note.
class DiracKinetic : public virtual rStatic_HT, private rStatic_HT_Imp
{
public:
    virtual void          GetEnergy(EnergyBreakdown&,const rDM_CD* cd ) const;
    virtual std::ostream& Write    (std::ostream&) const;
    virtual bool          IsRelativistic() const {return true;}
private:
    virtual rsmat_t MakeMatrix(const robs_t*,const Spin&) const;
};

class RestMass : public virtual rStatic_HT, private rStatic_HT_Imp
{
public:
    virtual void          GetEnergy(EnergyBreakdown&,const rDM_CD* cd ) const;
    virtual std::ostream& Write    (std::ostream&) const;
    virtual bool          IsRelativistic() const {return true;}
private:
    virtual rsmat_t MakeMatrix(const robs_t*,const Spin&) const;
};

// The ion-ion (nuclear-nuclear) repulsion energy is now the T-templated IonIon<T>
// (qchem.Hamiltonian.Internal.IonIon); the molecular Hamiltonians build IonIon<double>.

//
//  Electron-Nuclear attraction potential.
//
class Ven : public virtual rStatic_HT, private rStatic_HT_Imp
{
public:
    typedef std::shared_ptr<const Structure> st_t;
    Ven(const st_t& st);
    virtual void          GetEnergy(EnergyBreakdown&,const rDM_CD* cd) const;
    virtual std::ostream& Write    (std::ostream&) const;
private:
    virtual rsmat_t MakeMatrix(const robs_t*,const Spin&) const;
    st_t theStructure;
};

//###############################################################################
//
//  Local pseudopotential electron-ion term: the pseudized replacement for Ven.  Instead of the analytic
//  -Z/r nuclear attraction it quadratures the smooth real-space V_loc(r) on the molecular/atomic mesh,
//  <chi_i|V_loc|chi_j> = Sum_g w_g chi_i(r_g) chi_j(r_g) V_loc(r_g) (= the XC-path MatrixOverlap shape).
//  STATIC (density-independent), so it is built once.  V_loc is the species radial field's real-space view (SpeciesRadialField::ValueR).
//
class PP_Local : public virtual rStatic_HT, private rStatic_HT_Imp
{
public:
    //! The virial theorem needs a Coulombic (degree -1 homogeneous) potential; a pseudopotential is not
    //! (erf-screened local part + KB projectors), so the SCF drops both the virial gate and column (V1.27).
    virtual bool IsVirialValid() const {return false;}
    typedef std::shared_ptr<const Structure> st_t;
    typedef std::shared_ptr<const BasisSet::SpeciesRadialField> vloc_t;
    PP_Local(const st_t& st, vloc_t vloc, const qcMesh::MeshParams& mp);
    virtual void          GetEnergy(EnergyBreakdown&,const rDM_CD* cd) const;   // Een (PP local) = DM_Contract
    virtual std::ostream& Write    (std::ostream&) const;
private:
    virtual rsmat_t MakeMatrix(const robs_t*,const Spin&) const;
    st_t             theStructure;
    vloc_t           itsVloc;
    qcMesh::MeshParams itsMeshParams;
};

//###############################################################################
//
//  Separable (Kleinman-Bylander) NON-LOCAL pseudopotential term.  V_NL = Sum_{a,p} D_p Sum_m |beta_p Y_lm>
//  <beta_p Y_lm| -- real symmetric, STATIC.  This is the repulsive (for the occupied valence l-channels)
//  piece that lifts the over-bound local-only spectrum back to the all-electron valence eigenvalues.
//
//  TWO assemblies, chosen by a capability cross-cast on the block (the angular factor decides):
//   * EXPLICIT-angular (molecular/Cartesian) -- the 3-D mesh route: per atom, per projector, per m=-l..l a
//     rank-1 D|b><b| with b_i = <chi_i|beta_p(|r-R|) Y_lm(rhat)> (mesh quadrature).
//   * IMPLICIT-angular (ATOMIC, radial: BasisSet::ImplicitAngular_IBS) -- the per-l RADIAL route, because a
//     radial block's stored chi_i OMIT the irrep's Y_lm and the 3-D route would silently give nonsense
//     (an l=0 projector leaking into every l block; every l>=1 projector integrating to zero -- the
//     occupied-d defect found via MnO, doc/SymmetryUpgradePlan.md sec 7 step 7).  See MakeMatrixRadial.
//  Verified against the reciprocal (2l+1)P_l form; per-l gate A_PP.PerLKleinmanBylanderOracle (s,p,d,f).
//
class PP_NonLocal : public virtual rStatic_HT, private rStatic_HT_Imp
{
public:
    //! The virial theorem needs a Coulombic (degree -1 homogeneous) potential; a pseudopotential is not
    //! (erf-screened local part + KB projectors), so the SCF drops both the virial gate and column (V1.27).
    virtual bool IsVirialValid() const {return false;}
    typedef std::shared_ptr<const Structure> st_t;
    typedef std::shared_ptr<const BasisSet::SpeciesProjectorSet_R> sep_t;
    PP_NonLocal(const st_t& st, sep_t sep, const qcMesh::MeshParams& mp);
    virtual void          GetEnergy(EnergyBreakdown&,const rDM_CD* cd) const;   // Een (PP nonlocal) = DM_Contract
    virtual std::ostream& Write    (std::ostream&) const;
private:
    virtual rsmat_t MakeMatrix(const robs_t*,const Spin&) const;
    //! The ATOMIC (implicit-angular / radial) assembly -- see the body: a radial block's stored functions
    //! omit the irrep's Y_lm, so the 3-D mesh route cannot form <chi|beta Y_lm> and must not be used.
    rsmat_t MakeMatrixRadial(const BasisSet::ImplicitAngular_IBS&, size_t n) const;
    st_t             theStructure;
    sep_t            itsSep;
    qcMesh::MeshParams itsMeshParams;
};

//###############################################################################
//
//  Implementation of the Coulomb potential
//
//            /
// Vee(r_1) = | Ro(r_2)/r_12 d^3 r_2
//           /
//
// Ro is the exact charge density calculated from sum(Dab*a*b) using the density
// matrix and orbital basis functions.  This is the coulomb potential used in Hartree-Fock
// calculations.
//
// Shared whole-system machinery for the 4-index HF terms (Coulomb/exchange).  Owns the version guard, the
// composite-basis walk, the per-irrep block cache, and the Fock-build/energy plumbing that Vee and Vxc used
// to duplicate.  A concrete term supplies ONLY the one line that differs -- which canonical-pair contraction
// to run (AccumulateAll: Direct vs Exchange) -- plus an optional Fock Scale (Vxc's K coefficient).  Mirrors
// the tDynamic_HT / tDynamic_HT_Imp interface/impl split, so the tDynamic_HF_HT interface itself stays
// data-free.
class Dynamic_HF_HT_Imp : public virtual rDynamic_HF_HT
{
public:
    //! Fock build: assemble the whole-system blocks ONCE per (density, spin) from the composite \a wholeBasis
    //! using ERI4 bra-ket symmetry (canonical pairs -> ScatterBoth), cache the per-irrep blocks, return this
    //! irrep's block.  \a wholeBasis is required (HF is whole-system); a null basis throws.
    virtual const rsmat_t& GetMatrix(const robs_t*,const Spin&,const rChargeDensity*,const rbs_t* wholeBasis) const;
protected:
    //! The one operation that distinguishes Coulomb from exchange: scatter \a dm across canonical irrep pairs
    //! into the zeroed per-irrep blocks \a X (one per irrep, same order as the density's leaves).
    virtual void   AccumulateAll(std::vector<rsmat_t>& X,const rDM_CD* dm) const=0;
    //! \name THE SPIN AXIS (V1.37 step 3) -- three questions, answered by the concrete term, that used to
    //! be answered by WHICH TYPE it was (Vee / Vxc(-1/2) / two Vxc(-1) inside a VxcPol).
    //!@{
    //! Which spin keys my blocks: \c Spin::None when the operator ignores the channel (Coulomb sees the
    //! total), \a s itself when it is per channel (exchange is same-spin).
    virtual Spin           CacheSpin (const Spin& s) const=0;
    //! The density my spin-\a s blocks are built FROM: the total, or the \a s channel of \a cd.
    virtual const rDM_CD*  DensityFor(const rChargeDensity* cd, const Spin& s) const=0;
    //! Fock coefficient applied to every block after the scatter (1 for Coulomb; the K coefficient for Vxc,
    //! which depends on whether the block is a folded doublet or one channel).
    virtual double         Scale(const Spin& s) const=0;
    //!@}
    //! Contract \a cd into the whole-system blocks for spin \a s (a CacheSpin) if stale for this density.
    //! Uses itsWholeBasis (stashed from the Fock build), so GetEnergy -- which has no whole-basis -- gets the
    //! same symmetry-banked contraction for its (post-diagonalization) density.
    const std::map<std::string,rsmat_t>& ContractAll(const rChargeDensity* cd, const Spin& s) const;

    mutable const rbs_t* itsWholeBasis=nullptr;    //!< whole basis (stashed from the Fock build; stable across the run)
    //! The J (Coulomb) or K (exchange) per-irrep blocks of ONE spin: accumulated (over irreps) and contracted
    //! (over that spin's D) for the density ID'd by \c version.  Keyed by ab-basis BasisSetID, already scaled.
    struct Blocks { size_t version=size_t(-1); std::map<std::string,rsmat_t> jk; };
    mutable std::map<Spin,Blocks> itsJKs;          //!< one Blocks per CacheSpin this run asks for
};

class Vee : public Dynamic_HF_HT_Imp
{
public:
    virtual void          GetEnergy(EnergyBreakdown&,const rDM_CD* cd ) const;
    virtual std::ostream& Write    (std::ostream&) const;
protected:
    virtual void          AccumulateAll(std::vector<rsmat_t>& X,const rDM_CD* dm) const;
    // Coulomb sees the TOTAL density and ignores the channel: one set of blocks serves every spin.
    virtual Spin          CacheSpin (const Spin&) const {return Spin::None;}
    virtual const rDM_CD* DensityFor(const rChargeDensity* cd, const Spin&) const;
    virtual double        Scale(const Spin&) const {return 1.0;}
};

//###############################################################################
//
//  Hartree-Fock exchange -- ONE term for either imposed spin subgroup (V1.37 step 3).
//
//  Exchange is SAME-SPIN: \f$F^\sigma\mathrel{+}=-K[D_\sigma]\f$.  The block's spin says what \f$D_\sigma\f$ is:
//    * Up/Down (imposed U(1)_z): that channel of the density, coefficient -1;
//    * None (imposed SU(2), the folded doublet): the whole density IS \f$2D_\sigma\f$, so \f$-K[D_\sigma]=
//      -\tfrac12K[D_{tot}]\f$ -- the RHF \f$-\tfrac12\f$, read off the label instead of a second type.
//  (This replaced Vxc(-1/2) + a VxcPol forwarding to two Vxc(-1): the coefficient and the density were the
//  only differences, and both are functions of the block's spin.)
//
class Vxc : public Dynamic_HF_HT_Imp
{
public:
    virtual void           GetEnergy(EnergyBreakdown&,const rDM_CD* cd ) const;
    virtual std::ostream&  Write    (std::ostream&) const;
protected:
    virtual void          AccumulateAll(std::vector<rsmat_t>& X,const rDM_CD* dm) const;
    virtual Spin          CacheSpin (const Spin& s) const {return s;}          // same-spin: per channel
    virtual const rDM_CD* DensityFor(const rChargeDensity* cd, const Spin& s) const;
    virtual double        Scale(const Spin& s) const {return s==Spin::None ? -0.5 : -1.0;}
};

//###############################################################################
//
//  Implementation of the Coulomb potential
//
//            /
// Vee(r_1) = | Ro_fit(r_2)/r_12 d^3 r_2
//           /
//
// Where Ro is actually a fitted charge density.  This is the potential that is typically
// used in DFT calculations.  Ro_fit is expanded in a auxilliary basis set. The matrix elements
// involve three center integrals hence avoiding the four center integrals encountered in
// a Hartree-Fock calculation.
//
class FittedVee : public virtual rDynamic_HT, private rDynamic_HT_Imp
{
public:
    typedef std::shared_ptr<const BasisSet::rFIT_CD_ABS> fbs_t;   //!< the charge-density (Coulomb-metric) fit face
    FittedVee(fbs_t& chargeDensityFitBasisSet, double numElectrons);
    ~FittedVee();   // anchored in the Imp TU (FittedCD complete there) so the unique_ptr can delete it
    //! \copydoc tDynamic_HT::RefreshForDensity
    //! The fitted charge density is k-INDEPENDENT and refit once per density serial, so it is exactly what
    //! this phase is for.  It used to be fit LAZILY inside \c MakeMatrix, i.e. from inside the block loop.
    virtual void RefreshForDensity(const rChargeDensity* cd) const override;
    virtual void          GetEnergy(EnergyBreakdown&,const rDM_CD* cd) const;
    virtual std::ostream& Write    (std::ostream& os) const {return os;}
private:
    virtual rsmat_t MakeMatrix(const robs_t*,const Spin&,const rChargeDensity* cd) const;
    std::unique_ptr<ChargeDensity::FittedCD> itsFittedChargeDensity;   //!< owned (was a leaked raw ptr)
};


//###############################################################################
//
//  Linear least squares fit of the exchange-correlation potential -- ONE term for either imposed spin
//  subgroup (V1.37 step 3; it replaced FittedVxc + FittedVxcPol + FittedVcorrPol).  The fit basis set is
//  inserted by the constructor and is not owned; the XC functional IS owned (shared) here.
//
//  SPIN-NATIVE THROUGHOUT.  The functional is consumed through its two-channel face,
//  v^sigma(rho_up, rho_dn) and eps^sigma(rho_up, rho_dn), and the term supplies the channel densities:
//    * a polarized density hands over its Up/Down channels;
//    * a density that resolves no spin -- the folded doublet of an SU(2) run, or the spin-agnostic seed of a
//      polarized one -- hands over rho/2 for both, which is the exact zeta=0 collapse.
//  The imposed subgroup enters ONCE, at construction: it says which spin irreps this term will be asked for,
//  so a fitter pair exists for each of them before the first block loop (no lazy insertion, R1.0h).
//
//  TWO fits per spin irrep, on the SAME fit basis (so the 3-centre integrals are computed once):
//    V (GetMatrix / MakeMatrix) fits the POTENTIAL v^sigma        -> the Fock/KS block of that spin.
//    E (GetEMatrix)             fits the ENERGY DENSITY eps^sigma -> E_xc = Sum_sigma Tr(D_sigma <i|eps^sigma|j>).
//  They are genuinely different matrices (v = eps + rho d(eps)/d(rho); a factor 4/3 for Slater exchange),
//  which is the whole reason tDynamic_CC's energy face is named GetEMatrix rather than GetMatrix (V1.3).
//  eps^sigma is PER CHANNEL because exchange's energy density is (eps_x(rho_up) differs from eps_x(rho_dn));
//  correlation's is the same for both, and a composite sums them (ExFunctional::GetEpsXc(up,dn,s)).
//
//  The v fit keys on the Fock pass's density (rho_in), the eps fit on the energy pass's (rho_out), so the
//  two guards cannot be shared -- and neither can be a lazy insert: both slots exist per spin from
//  construction and are only ever REFILLED.
//
class FittedVxc : public virtual rDynamic_HT, private rDynamic_HT_Imp
{
public:
    typedef std::shared_ptr<const BasisSet::rFIT_SF_ABS> fbs_t;   //!< the scalar-function (overlap-metric) fit face
    typedef std::shared_ptr<      ExFunctional>  ex_t;

    //! \a g names the imposed spin subgroup: the spin irreps this term will be asked for (None; or Up+Down).
    FittedVxc(fbs_t&, ex_t&, SpinGroup g);
    ~FittedVxc();
    //! \copydoc tDynamic_HT::RefreshForDensity
    //! THE EAGER PHASE (R1.0h): the v^sigma fits are k-INDEPENDENT, so hoist them out of the block loop --
    //! one per spin irrep, each from BOTH channel densities of \a cd.  ⚠ Only the V half: the eps fits key
    //! on the ENERGY pass's density, not this pass's (see GetEMatrix), so warming them here would fit the
    //! wrong one and they would be refit anyway.
    virtual void RefreshForDensity(const rChargeDensity* cd) const override;
    virtual void          GetEnergy       (EnergyBreakdown&,const rDM_CD*) const override;
    //! The ENERGY block of spin \a s: re-fits eps^s for this density and returns Sum_a c_a <Oi|f_a|Oj>.
    virtual const rsmat_t& GetEMatrix(const robs_t*,const Spin&,const rChargeDensity* cd) const override;
    virtual std::ostream& Write           (std::ostream&) const override;
private:
    virtual rsmat_t MakeMatrix(const robs_t*,const Spin&,const rChargeDensity*) const override;

    //! The fitter PAIR of one spin irrep, each half guarded by the density serial it currently holds.
    struct SpinFit
    {
        std::unique_ptr<Fitting::FunctionFitter_Scalar> v;     //!< V: the v^sigma fit (the Fock pass's density)
        std::unique_ptr<Fitting::FunctionFitter_Scalar> eps;   //!< E: the eps^sigma fit (the energy pass's density)
        size_t  vVersion  =size_t(-1);
        size_t  epsVersion=size_t(-1);
        rsmat_t epsMat;                                        //!< GetEMatrix's returned block (per spin: no aliasing)
    };
    //! This spin irrep's pair -- THROWS for a spin the term was not built for (a block of the other subgroup).
    SpinFit& FitFor(const Spin& s) const;
    //! Bring \a s's v fit up to \a cd (a no-op when it already holds this serial).
    void EnsureVFit(SpinFit&, const Spin& s, const rChargeDensity* cd) const;

    ex_t                                itsEx;     //!< the XC functional (owned, shared): both faces, both fits
    SpinGroup                           itsGroup;  //!< the imposed subgroup this term was built for
    mutable std::map<Spin,SpinFit>      itsFits;   //!< one pair per spin irrep of the subgroup, from construction
};

} //namespace