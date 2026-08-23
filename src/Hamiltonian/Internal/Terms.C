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
import qchem.Pseudopotential.LocalPotential;   // LocalPotential_R (the real-space PP local view)
import qchem.Pseudopotential.SeparablePotential; // SeparablePotential_R (the real-space KB projector view)
import qchem.Mesh;                             // qcMesh::MeshParams (the quadrature mesh spec)
import qchem.BasisSet.ImplicitAngular_IBS;      // the radial/implicit-Y_lm capability (atomic KB route)


export namespace qchem::Hamiltonian
{

using ChargeDensity::Polarized_CD;

// The non-relativistic kinetic ENERGY term is now the T-templated Kinetic<T>
// (qchem.Hamiltonian.Internal.Kinetic); the molecular Hamiltonians build Kinetic<double>.

// Relativistic kinetic ENERGY term (Dirac \f$c\,\vec\sigma\cdot\vec p\f$); consumes the RKB-assembled
// relativistic kinetic block directly (no 1/2). See Imp/DiracKinetic.C for the unverified-factor note.
class DiracKinetic : public virtual rStatic_HT, private rStatic_HT_Imp
{
public:
    virtual void          GetEnergy(EnergyBreakdown&,const rDM_CD* cd ) const;
    virtual std::ostream& Write    (std::ostream&) const;
    virtual bool          IsPolarized   () const {return true;}
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
//  <chi_i|V_loc|chi_j> = Sum_g w_g chi_i(r_g) chi_j(r_g) V_loc(r_g) (= the XC-path WeightedOverlap shape).
//  STATIC (density-independent), so it is built once.  V_loc is the real-space PP face (LocalPotential_R).
//
class PP_Local : public virtual rStatic_HT, private rStatic_HT_Imp
{
public:
    //! The virial theorem needs a Coulombic (degree -1 homogeneous) potential; a pseudopotential is not
    //! (erf-screened local part + KB projectors), so the SCF drops both the virial gate and column (V1.27).
    virtual bool IsVirialValid() const {return false;}
    typedef std::shared_ptr<const Structure> st_t;
    typedef std::shared_ptr<const Pseudopotential::LocalPotential_R> vloc_t;
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
    typedef std::shared_ptr<const Pseudopotential::SeparablePotential_R> sep_t;
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
    //! Fock build: assemble the whole-system blocks ONCE per density from the composite \a wholeBasis using
    //! ERI4 bra-ket symmetry (canonical pairs -> ScatterBoth), cache the per-irrep blocks, return this
    //! irrep's block.  \a wholeBasis is required (HF is whole-system); a null basis throws.
    virtual const rsmat_t& GetMatrix(const robs_t*,const Spin&,const rChargeDensity*,const rbs_t* wholeBasis) const;
protected:
    //! The one operation that distinguishes Coulomb from exchange: scatter \a dm across canonical irrep pairs
    //! into the zeroed per-irrep blocks \a X (one per irrep, same order as the density's leaves).
    virtual void   AccumulateAll(std::vector<rsmat_t>& X,const rDM_CD* dm) const=0;
    //! Fock coefficient applied to every block after the scatter (1 for Coulomb; the K coefficient for Vxc).
    virtual double Scale() const {return 1.0;}
    //! Contract the density into the whole-system blocks itsJKs (keyed by BasisSetID) if stale for this
    //! density.  Uses itsWholeBasis (stashed from the Fock build), so GetEnergy -- which has no whole-basis --
    //! gets the same symmetry-banked contraction for its (post-diagonalization) density.
    void ContractAll(const rChargeDensity* cd) const;

    mutable size_t itsCD_Version=size_t(-1);      //!< ID number for the most recent charge density (CD)
    mutable const rbs_t* itsWholeBasis=nullptr;    //!< whole basis (stashed from the Fock build; stable across the run)
    //! The J (Coulomb) or K (exchange) per-irrep blocks: accumulated (over irreps) and contracted (over Dcd)
    //! for the current charge density ID'd by itsCD_Version.  Keyed by ab-basis BasisSetID, already scaled.
    mutable std::map<std::string,rsmat_t> itsJKs;
};

class Vee : public Dynamic_HF_HT_Imp
{
public:
    virtual void          GetEnergy(EnergyBreakdown&,const rDM_CD* cd ) const;
    virtual std::ostream& Write    (std::ostream&) const;
protected:
    virtual void AccumulateAll(std::vector<rsmat_t>& X,const rDM_CD* dm) const;
};

//###############################################################################
//
//  Hartree-Fock unpolarized and polarized exchange potentials.
//
class Vxc : public Dynamic_HF_HT_Imp
{
public:
    //! \a exchangeScale is the K coefficient in the Fock: -1/2 for the (spin-summed) RHF term, -1 for each
    //! spin channel of the polarized term (VxcPol owns two Vxc(-1)).  Explicit -- no hidden convention: the
    //! block contracts whatever density it is handed (total for RHF; a single spin channel for VxcPol).
    explicit Vxc(double exchangeScale) : itsScale(exchangeScale) {}
    virtual void           GetEnergy(EnergyBreakdown&,const rDM_CD* cd ) const;
    virtual std::ostream&  Write    (std::ostream&) const;
protected:
    virtual void   AccumulateAll(std::vector<rsmat_t>& X,const rDM_CD* dm) const;
    virtual double Scale() const {return itsScale;}   //!< the Fock K coefficient, applied to every block
private:
    const double itsScale;                          //!< K coefficient in the Fock (-1/2 RHF, -1 per-spin)
};

// Polarized HF exchange = two spin-channel Vxc(-1): dispatch per spin, feeding each its own spin density
// (K^sigma from D^sigma).  Mirrors FittedVxcPol's owned-pair structure -- keeps the fitted and HF polarized
// terms consistent.
class VxcPol : public virtual rDynamic_HF_HT
{
public:
    VxcPol();
   ~VxcPol();
    virtual void           GetEnergy(EnergyBreakdown&,const rDM_CD* cd ) const;
    virtual bool           IsPolarized() const {return true;}
    virtual std::ostream&  Write    (std::ostream&) const;
    virtual const rsmat_t& GetMatrix(const robs_t*,const Spin&,const rChargeDensity*,const rbs_t* wholeBasis) const;
private:
    Vxc* itsUpVxc  ;   //!< owned; spin-up exchange   (K coefficient -1)
    Vxc* itsDownVxc;   //!< owned; spin-down exchange (K coefficient -1)
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
    virtual void          GetEnergy(EnergyBreakdown&,const rDM_CD* cd) const;
    virtual std::ostream& Write    (std::ostream& os) const {return os;}
private:
    virtual rsmat_t MakeMatrix(const robs_t*,const Spin&,const rChargeDensity* cd) const;
    std::unique_ptr<ChargeDensity::FittedCD> itsFittedChargeDensity;   //!< owned (was a leaked raw ptr)
};


//###############################################################################
//
//  Linear least squares fit the unpolarized and polarized exchange-correlation potential.  The fit basis set
//  is inserted by the constructor and is not owned by FittedVxc; the XC functional IS owned (shared) here.
//
//  TWO fits, on the SAME fit basis (so the 3-centre integrals are computed once):
//    V (GetMatrix / MakeMatrix) fits the POTENTIAL v_xc(rho(r))       -> the Fock/KS block.
//    E (GetEMatrix)             fits the ENERGY DENSITY eps_xc(rho(r)) -> E_xc = integral eps_xc rho.
//  They are genuinely different matrices (v_xc = eps_xc + rho d(eps_xc)/d(rho); a factor 4/3 for Slater
//  exchange), which is the whole reason tDynamic_CC's energy face is named GetEMatrix rather than GetMatrix
//  -- before that, delivering the second matrix needed a SEPARATE rDynamic_CC object (the retired
//  FittedEpsXc adapter) hung off this term purely to carry it (V1.3).
//  Uniform for exchange AND correlation: the E fit reads the functional's own eps_xc (eps_x = 3/4 v_x for
//  Dirac exchange; eps_c != 3/4 v_c for correlation; eps_xc from libxc), so no functional needs the 3/4
//  special case -- this retires the old 3/4-virial exchange shortcut (which broke for gradient functionals).
//
class FittedVxc : public virtual rDynamic_HT, private rDynamic_HT_Imp
{
public:
    typedef std::shared_ptr<const BasisSet::rFIT_SF_ABS> fbs_t;   //!< the scalar-function (overlap-metric) fit face
    typedef std::shared_ptr<ExFunctional>     ex_t;

    FittedVxc(fbs_t& VxcFitBasisSet, ex_t&);
    ~FittedVxc();
    virtual void          GetEnergy       (EnergyBreakdown&,const rDM_CD*) const override;
    //! The ENERGY block: re-fits eps_xc for this density and returns Sum_a c_a <Oi|f_a|Oj> for contraction.
    virtual const rsmat_t& GetEMatrix(const robs_t*,const Spin&,const rChargeDensity* cd) const override;
    virtual std::ostream& Write           (std::ostream&) const override;
private:
    virtual rsmat_t MakeMatrix(const robs_t*,const Spin&,const rChargeDensity*) const override;

    ex_t itsEx;   //!< the XC functional (owned, shared): supplies BOTH GetVxc (V) and GetEpsXc (E)
    std::unique_ptr<Fitting::FunctionFitter_Scalar> itsFitter;    //!< V: the v_xc fit
    std::unique_ptr<Fitting::FunctionFitter_Scalar> itsEpsFitter; //!< E: the eps_xc fit (same fit basis)
    mutable rsmat_t     itsEpsMat;                   //!< GetEMatrix's returned block
    mutable size_t      itsEpsVersion=size_t(-1);    //!< density serial the eps_xc fit was last computed for
};

//! Polarized exchange as a pure FORWARDER to its two per-spin children -- the same shape as \c VxcPol,
//! its Hartree-Fock twin (R2.19).  It owns no matrix of its own: each child is a full caching term, so
//! this class hands back the CHILD'S cached reference rather than copying it into scratch.  That is why
//! it does NOT derive from \c rDynamic_HT_Imp_NoCache: it has nothing to compute, so it has no
//! \c MakeMatrix, and a scratch slot existed only to have something to return a reference to.
class FittedVxcPol : public virtual rDynamic_HT
{
public:
    typedef std::shared_ptr<const BasisSet::rFIT_SF_ABS> fbs_t;   //!< the scalar-function (overlap-metric) fit face
    typedef std::shared_ptr<      ExFunctional>  ex_t;

    FittedVxcPol(fbs_t&, ex_t&);
   ~FittedVxcPol();
    //! Forward to this spin's child and return ITS cached block.  Reference lifetime is therefore the
    //! child's: valid until the next call for the same Irrep on that child.
    virtual const rsmat_t& GetMatrix(const robs_t*,const Spin&,const rChargeDensity* cd) const;
    // Required by HamiltonianTerm
    virtual void GetEnergy       (EnergyBreakdown&,const rDM_CD* cd         ) const;
    virtual bool IsPolarized() const {return true;}

    virtual std::ostream&   Write(std::ostream&) const;
private:
    rDynamic_HT* itsUpVxc  ; //Spin up.
    rDynamic_HT* itsDownVxc; //Spin down.

};

//###############################################################################
//
//  Polarized (spin-native) correlation term.  Unlike FittedVxcPol -- which delegates to two INDEPENDENT
//  single-channel FittedVxc, valid only because Slater exchange is channel-separable -- correlation
//  v_c^sigma(rho_up,rho_down) COUPLES both channels (through r_s and zeta), so this term fits the
//  SpinCorrelation functional against the FULL Polarized_CD at each mesh point.  The Fock build calls
//  MakeMatrix per spin (each fits v_c^sigma); the energy E_c = integral eps_c(rho_up,rho_down) rho uses a
//  SECOND eps_c fit on the same fit basis (GetEMatrix -- the E face of the V/E pair) that the polarized
//  density contracts over both channels.  (That fit used to live in a separate rDynamic_CC adapter,
//  FittedEpsCPol -- a clone of FittedVxc's FittedEpsXc; both died with the GetEMatrix split, V1.3.)  The seed
//  iteration (a spin-agnostic total density, not yet a Polarized_CD) collapses to v_c^P(rho) via
//  rho_up=rho_down=rho/2 -- the same robustness FittedVxcPol needed (cd85d13c).
//
class FittedVcorrPol : public virtual rDynamic_HT, private rDynamic_HT_Imp_NoCache
{
public:
    typedef std::shared_ptr<const BasisSet::rFIT_SF_ABS> fbs_t;   //!< the scalar-function (overlap-metric) fit face
    typedef std::shared_ptr<SpinCorrelation>            corr_t;  //!< the spin-native correlation functional

    FittedVcorrPol(fbs_t&, corr_t&);
   ~FittedVcorrPol();
    virtual void GetEnergy (EnergyBreakdown&, const rDM_CD* cd) const override;
    //! The ENERGY block: fits eps_c(rho_up,rho_down) from the full Polarized_CD and returns its overlap
    //! matrix.  Spin-INDEPENDENT as a value, so contracting it over both channels gives
    //! E_c = integral eps_c (rho_up+rho_down).
    virtual const rsmat_t& GetEMatrix(const robs_t*, const Spin&, const rChargeDensity* cd) const override;
    virtual bool IsPolarized() const override {return true;}
    virtual std::ostream& Write(std::ostream&) const override;
private:
    virtual rsmat_t MakeMatrix(const robs_t*, const Spin&, const rChargeDensity* cd) const override;

    corr_t itsCorr;                                                      //!< the correlation functional (owned)
    std::unique_ptr<Fitting::FunctionFitter_Scalar> itsVcFitter; //!< v_c^sigma potential fit
    std::unique_ptr<Fitting::FunctionFitter_Scalar> itsEpsFitter;//!< E_c = integral eps_c rho (energy)
    mutable rsmat_t                                         itsEpsMat;   //!< GetEMatrix's returned block
};


} //namespace