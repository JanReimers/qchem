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
import qchem.Hamiltonian.Internal.LDAVxc;     // LDAVxc (the v_xc fit client held by FittedVxc)
import qchem.Structure;
import qchem.Fitting.FunctionFitter;          // Fitting::FunctionFitter (composed; clients never see the impl)
import qchem.ChargeDensity;
import qchem.FittedCD;
import qchem.Hamiltonian.Types;
import qchem.Pseudopotential.LocalPotential;   // LocalPotential_R (the real-space PP local view)
import qchem.Pseudopotential.SeparablePotential; // SeparablePotential_R (the real-space KB projector view)
import qchem.Mesh;                             // qcMesh::MeshParams (the quadrature mesh spec)


export namespace qchem::Hamiltonian
{

using ChargeDensity::Polarized_CD;

// Non-relativistic kinetic ENERGY term \f$T=-\tfrac12\nabla^2 = \tfrac12\langle p^2\rangle\f$.
// (Applies the 1/2 to the basis-set \f$\langle p^2\rangle\f$ block; see Imp/Kinetic.C.)
class Kinetic : public virtual rStatic_HT, private rStatic_HT_Imp
{
public:
    virtual void          GetEnergy(EnergyBreakdown&,const rDM_CD* cd ) const;
    virtual std::ostream& Write    (std::ostream&) const;
private:
    virtual rsmat_t CalculateMatrix(const robs_t*,const Spin&) const;
};

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
    virtual rsmat_t CalculateMatrix(const robs_t*,const Spin&) const;
};

class RestMass : public virtual rStatic_HT, private rStatic_HT_Imp
{
public:
    virtual void          GetEnergy(EnergyBreakdown&,const rDM_CD* cd ) const;
    virtual std::ostream& Write    (std::ostream&) const;
    virtual bool          IsRelativistic() const {return true;}
private:
    virtual rsmat_t CalculateMatrix(const robs_t*,const Spin&) const;
};

//
//  Nuclear-Nuclear repulsion energy.  It does not *need* to contribute to the Hamiltonian matrix (a
//  constant), only the Enn energy: a direct Coulomb pair sum for a finite molecule, an Ewald lattice sum
//  for a periodic cell (NuclearRepulsion picks by Structure::isFinite()).  The ion charge of each atom is
//  supplied by a Z->charge callback: the all-electron default is identity (itsZ), a pseudopotential maps
//  Z -> its valence Zion core charge -- so the SAME term serves both, mirroring the PW PW_IonIon term.
//
class Vnn : public virtual rStatic_HT, private rStatic_HT_Imp
{
public:
    typedef std::shared_ptr<const Structure> st_t;
    Vnn(const st_t& st);                                            //!< all-electron: ion charge = itsZ
    Vnn(const st_t& st, std::function<double(int)> zionOf);        //!< pseudopotential: Z -> Zion core charge
    virtual void          GetEnergy(EnergyBreakdown&,const rDM_CD* cd) const;
    virtual std::ostream& Write    (std::ostream&) const;
private:
    virtual rsmat_t CalculateMatrix(const robs_t*,const Spin&) const;
    st_t theStructure;
    std::function<double(int)> itsZionOf;   //!< Z -> ion core charge (identity for the all-electron baseline)
};

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
    virtual rsmat_t CalculateMatrix(const robs_t*,const Spin&) const;
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
    typedef std::shared_ptr<const Structure> st_t;
    typedef std::shared_ptr<const Pseudopotential::LocalPotential_R> vloc_t;
    PP_Local(const st_t& st, vloc_t vloc, const qcMesh::MeshParams& mp);
    virtual void          GetEnergy(EnergyBreakdown&,const rDM_CD* cd) const;   // Een (PP local) = DM_Contract
    virtual std::ostream& Write    (std::ostream&) const;
private:
    virtual rsmat_t CalculateMatrix(const robs_t*,const Spin&) const;
    st_t             theStructure;
    vloc_t           itsVloc;
    qcMesh::MeshParams itsMeshParams;
};

//###############################################################################
//
//  Separable (Kleinman-Bylander) NON-LOCAL pseudopotential term.  Per atom, per projector p (angular
//  momentum l, strength D = Coefficient), per m=-l..l it is the rank-1 outer product D|b><b| with the
//  projection vector  b_i = <chi_i | beta_p(|r-R|) Y_lm(rhat)>  (mesh quadrature; beta_p = BetaR, Y_lm
//  unit-normalised on the sphere).  V_NL = Sum_{a,p,m} D |b><b| -- real symmetric, STATIC.  This is the
//  repulsive (for the occupied valence l-channels) piece that lifts the over-bound local-only spectrum
//  back to the all-electron valence eigenvalues.  Verified against the reciprocal (2l+1)P_l form.
//
class PP_NonLocal : public virtual rStatic_HT, private rStatic_HT_Imp
{
public:
    typedef std::shared_ptr<const Structure> st_t;
    typedef std::shared_ptr<const Pseudopotential::SeparablePotential_R> sep_t;
    PP_NonLocal(const st_t& st, sep_t sep, const qcMesh::MeshParams& mp);
    virtual void          GetEnergy(EnergyBreakdown&,const rDM_CD* cd) const;   // Een (PP nonlocal) = DM_Contract
    virtual std::ostream& Write    (std::ostream&) const;
private:
    virtual rsmat_t CalculateMatrix(const robs_t*,const Spin&) const;
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
    virtual rsmat_t CalcMatrix(const robs_t*,const Spin&,const rChargeDensity* cd) const;
    std::unique_ptr<ChargeDensity::FittedCD> itsFittedChargeDensity;   //!< owned (was a leaked raw ptr)
};


//###############################################################################
//
//  Linear least squares fit the unpolarized and polarized exchange-correlation potential.  The fit basis set is inserted by the
//  constructor and is not owned by FittedVxc.  The XC functional is owned by the inner LDAVxc.
//  Energy uses the exchange virial E_xc = 3/4 <rho|Vxc> -- correct for (pure) exchange; for correlation
//  use FittedVcorr below.
//
class FittedVxc : public virtual rDynamic_HT, private rDynamic_HT_Imp
{
public:
    typedef std::shared_ptr<const BasisSet::FIT_SF_ABS> fbs_t;   //!< the scalar-function (overlap-metric) fit face
    typedef std::shared_ptr<ExFunctional>     ex_t;

    FittedVxc(fbs_t& VxcFitBasisSet, ex_t&);
    ~FittedVxc();
    virtual void          GetEnergy       (EnergyBreakdown&,const rDM_CD*) const;
    virtual void          UseChargeDensity(const rChargeDensity*);
    virtual std::ostream& Write           (std::ostream&) const;
private:
    virtual rsmat_t CalcMatrix(const robs_t*,const Spin&,const rChargeDensity*) const;

    std::unique_ptr<Fitting::FunctionFitter_Scalar<double>> itsFitter; //!< COMPOSED v_xc fit (was inherited)
    LDAVxc* itsLDAVxc;   //!< the v_xc=Vxc(rho) function to fit (concrete: it IS the ScalarFFClient)
};

class FittedVxcPol : public virtual rDynamic_HT, private rDynamic_HT_Imp_NoCache
{
public:
    typedef std::shared_ptr<const BasisSet::FIT_SF_ABS> fbs_t;   //!< the scalar-function (overlap-metric) fit face
    typedef std::shared_ptr<      ExFunctional>  ex_t;

    FittedVxcPol(fbs_t&, ex_t&);
   ~FittedVxcPol();
    // Required by HamiltonianTerm
    virtual void GetEnergy       (EnergyBreakdown&,const rDM_CD* cd         ) const;
    virtual bool IsPolarized() const {return true;}

    virtual std::ostream&   Write(std::ostream&) const;
private:
    virtual rsmat_t CalcMatrix(const robs_t*,const Spin&,const rChargeDensity* cd) const;

    rDynamic_HT* itsUpVxc  ; //Spin up.
    rDynamic_HT* itsDownVxc; //Spin down.

};

//###############################################################################
//
//  This is a total energy term, not a matrix values Hamiltonian term.
//  A dedicated least-squares fit of an XC ENERGY DENSITY eps_xc(rho(r)) to the auxiliary fit basis,
//  exposed as a rDynamic_CC so the density contracts it:  E_xc = integral eps_xc rho = <rho|eps_xc_fit>.
//  Separate from a potential fit (different coefficients) but meant to SHARE its fit basis, so the
//  3-centre integrals are computed once.  Needed because eps_c != 3/4 v_c -- the exchange virial
//  (eps_x = 3/4 v_x) used by FittedVxc::GetEnergy is correct for exchange but wrong for correlation.
//  Composes a Fitting::FunctionFitter (from the Factory); Overlap is queried on it.
//
class FittedEpsXc : public virtual ChargeDensity::rDynamic_CC
{
public:
    typedef std::shared_ptr<const BasisSet::FIT_SF_ABS> fbs_t;   //!< the scalar-function (overlap-metric) fit face

    FittedEpsXc(fbs_t& fitBasisSet, const ExFunctional* ex);
    //! Re-fits eps_xc for this density and returns its matrix Sum_a c_a <Oi|f_a|Oj> for contraction.
    virtual const rsmat_t& GetMatrix(const robs_t*,const Spin&,const rChargeDensity* cd) const;
private:
    std::unique_ptr<Fitting::FunctionFitter_Scalar<double>> itsFitter;  //!< COMPOSED fitter (not inherited)
    const ExFunctional* itsEx;   //!< non-owning; the XC functional supplying eps_xc (owned by the term)
    mutable rsmat_t     itsMat;
};
//###############################################################################
//
//  Correlation term.  Same potential->matrix machinery as FittedVxc (fits v_c into H), but the ENERGY is
//  E_c = integral eps_c rho via a dedicated FittedEpsXc on the SAME fit basis -- NOT the 3/4 exchange
//  virial (eps_c != 3/4 v_c).  Pair with an exchange FittedVxc(Dirac) to assemble a full LDA Hamiltonian.
//
class FittedVcorr : public FittedVxc
{
public:
    FittedVcorr(fbs_t& VcorrFitBasisSet, ex_t&);
    virtual void GetEnergy(EnergyBreakdown&,const rDM_CD* cd) const;
private:
    FittedEpsXc itsEpsC;   //!< dedicated eps_c fit for the correlation energy (shares the fit basis)
};

class FittedEpsCPol;   // the polarized eps_c contraction client (defined in Imp/FittedVcorrPol.C)

//###############################################################################
//
//  Polarized (spin-native) correlation term.  Unlike FittedVxcPol -- which delegates to two INDEPENDENT
//  single-channel FittedVxc, valid only because Slater exchange is channel-separable -- correlation
//  v_c^sigma(rho_up,rho_down) COUPLES both channels (through r_s and zeta), so this term fits the
//  SpinCorrelation functional against the FULL Polarized_CD at each mesh point.  The Fock build calls
//  CalcMatrix per spin (each fits v_c^sigma); the energy E_c = integral eps_c(rho_up,rho_down) rho uses a
//  dedicated eps_c fit (FittedEpsCPol) that the polarized density contracts over both channels.  The seed
//  iteration (a spin-agnostic total density, not yet a Polarized_CD) collapses to v_c^P(rho) via
//  rho_up=rho_down=rho/2 -- the same robustness FittedVxcPol needed (cd85d13c).
//
class FittedVcorrPol : public virtual rDynamic_HT, private rDynamic_HT_Imp_NoCache
{
public:
    typedef std::shared_ptr<const BasisSet::FIT_SF_ABS> fbs_t;   //!< the scalar-function (overlap-metric) fit face
    typedef std::shared_ptr<SpinCorrelation>            corr_t;  //!< the spin-native correlation functional

    FittedVcorrPol(fbs_t&, corr_t&);
   ~FittedVcorrPol();
    virtual void GetEnergy (EnergyBreakdown&, const rDM_CD* cd) const;
    virtual bool IsPolarized() const {return true;}
    virtual std::ostream& Write(std::ostream&) const;
private:
    virtual rsmat_t CalcMatrix(const robs_t*, const Spin&, const rChargeDensity* cd) const;

    corr_t itsCorr;                                                      //!< the correlation functional (owned)
    std::unique_ptr<Fitting::FunctionFitter_Scalar<double>> itsVcFitter; //!< v_c^sigma potential fit
    std::unique_ptr<FittedEpsCPol>                          itsEpsC;     //!< E_c = integral eps_c rho (energy)
};


} //namespace