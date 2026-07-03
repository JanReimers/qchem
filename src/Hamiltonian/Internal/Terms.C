// Hamiltonian/Internal/Terms.C  Declare and export all Hamiltonian term types.
module;
#include <iosfwd>
#include <memory>
#include <map>
#include <string>
#include <vector>
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
class Kinetic : public virtual Static_HT, private Static_HT_Imp
{
public:
    virtual void          GetEnergy(EnergyBreakdown&,const DM_CD* cd ) const;
    virtual std::ostream& Write    (std::ostream&) const;
private:
    virtual rsmat_t CalculateMatrix(const obs_t*,const Spin&) const;
};

// Relativistic kinetic ENERGY term (Dirac \f$c\,\vec\sigma\cdot\vec p\f$); consumes the RKB-assembled
// relativistic kinetic block directly (no 1/2). See Imp/DiracKinetic.C for the unverified-factor note.
class DiracKinetic : public virtual Static_HT, private Static_HT_Imp
{
public:
    virtual void          GetEnergy(EnergyBreakdown&,const DM_CD* cd ) const;
    virtual std::ostream& Write    (std::ostream&) const;
    virtual bool          IsPolarized   () const {return true;}
    virtual bool          IsRelativistic() const {return true;}
private:
    virtual rsmat_t CalculateMatrix(const obs_t*,const Spin&) const;
};

class RestMass : public virtual Static_HT, private Static_HT_Imp
{
public:
    virtual void          GetEnergy(EnergyBreakdown&,const DM_CD* cd ) const;
    virtual std::ostream& Write    (std::ostream&) const;
    virtual bool          IsRelativistic() const {return true;}
private:
    virtual rsmat_t CalculateMatrix(const obs_t*,const Spin&) const;
};

//
//  Nuclear-Nuclear repulsion potential.  Actually this is just an energy, it does not *need* to make
//  a contribution to the Hamiltonian matrix.  For solids (structure type = Lattice_{1,2,3}D) this
//  should decide at run time to do an Ewald sum.  But and Ewald sum requires more information: ionic charges
//  around each nucleus.
//
class Vnn : public virtual Static_HT, private Static_HT_Imp
{
public:
    typedef std::shared_ptr<const Structure> st_t;
    Vnn(const st_t& st);
    virtual void          GetEnergy(EnergyBreakdown&,const DM_CD* cd) const;
    virtual std::ostream& Write    (std::ostream&) const;
private:
    virtual rsmat_t CalculateMatrix(const obs_t*,const Spin&) const;
    st_t theStructure;
};

//
//  Electron-Nuclear attraction potential.
//
class Ven : public virtual Static_HT, private Static_HT_Imp
{
public:
    typedef std::shared_ptr<const Structure> st_t;
    Ven(const st_t& st);
    virtual void          GetEnergy(EnergyBreakdown&,const DM_CD* cd) const;
    virtual std::ostream& Write    (std::ostream&) const;
private:
    virtual rsmat_t CalculateMatrix(const obs_t*,const Spin&) const;
    st_t theStructure;
};

//###############################################################################
//
//  Local pseudopotential electron-ion term: the pseudized replacement for Ven.  Instead of the analytic
//  -Z/r nuclear attraction it quadratures the smooth real-space V_loc(r) on the molecular/atomic mesh,
//  <chi_i|V_loc|chi_j> = Sum_g w_g chi_i(r_g) chi_j(r_g) V_loc(r_g) (= the XC-path WeightedOverlap shape).
//  STATIC (density-independent), so it is built once.  V_loc is the real-space PP face (LocalPotential_R).
//
class PP_Local : public virtual Static_HT, private Static_HT_Imp
{
public:
    typedef std::shared_ptr<const Structure> st_t;
    typedef std::shared_ptr<const Pseudopotential::LocalPotential_R> vloc_t;
    PP_Local(const st_t& st, vloc_t vloc, const qcMesh::MeshParams& mp);
    virtual void          GetEnergy(EnergyBreakdown&,const DM_CD* cd) const;   // Een (PP local) = DM_Contract
    virtual std::ostream& Write    (std::ostream&) const;
private:
    virtual rsmat_t CalculateMatrix(const obs_t*,const Spin&) const;
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
class PP_NonLocal : public virtual Static_HT, private Static_HT_Imp
{
public:
    typedef std::shared_ptr<const Structure> st_t;
    typedef std::shared_ptr<const Pseudopotential::SeparablePotential_R> sep_t;
    PP_NonLocal(const st_t& st, sep_t sep, const qcMesh::MeshParams& mp);
    virtual void          GetEnergy(EnergyBreakdown&,const DM_CD* cd) const;   // Een (PP nonlocal) = DM_Contract
    virtual std::ostream& Write    (std::ostream&) const;
private:
    virtual rsmat_t CalculateMatrix(const obs_t*,const Spin&) const;
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
class Vee : public virtual Dynamic_HT, private Dynamic_HT_Imp
{
public:
    virtual void          GetEnergy(EnergyBreakdown&,const DM_CD* cd ) const;
    virtual std::ostream& Write    (std::ostream&) const;
    //! Fock build WITH the whole (composite) basis: assemble the whole-system Coulomb ONCE per density using
    //! ERI4 bra-ket symmetry (canonical pairs -> ScatterBoth), cache the per-irrep blocks, and return this
    //! irrep's block.  With a null whole-basis (stand-alone tests) it falls back to the per-irrep base path.
    virtual const rsmat_t& GetMatrix(const obs_t*,const Spin&,const rChargeDensity*,const bs_t* wholeBasis) const;
    //! No-basis callers (the energy DM_Contract): reuse the whole-system cache when it is fresh for this
    //! density, else the per-irrep base path.
    virtual const rsmat_t& GetMatrix(const obs_t*,const Spin&,const rChargeDensity*) const;
private:
    virtual rsmat_t CalcMatrix(const obs_t*,const Spin&,const rChargeDensity* cd) const;
    //! (Re)build the whole-system Coulomb into itsJ (keyed by BasisSetID) if stale for this density.  Uses
    //! itsWholeBasis (stashed from the Fock build), so the energy path -- which has no whole-basis -- also
    //! gets the symmetry-banked build for its (post-diagonalization) density.
    void EnsureWholeSystem(const rChargeDensity* cd) const;
    mutable size_t itsAllVersion=size_t(-1);        //!< density serial the whole-system Coulomb was built for
    mutable const bs_t* itsWholeBasis=nullptr;      //!< whole basis (stashed from the Fock build; stable across the run)
    mutable std::map<std::string,rsmat_t> itsJ;     //!< per-irrep Coulomb blocks, keyed by ab-basis BasisSetID
};

//###############################################################################
//
//  Hartree-Fock unpolarized and polarized exchange potentials.
//
class Vxc : public virtual Dynamic_HT, private Dynamic_HT_Imp
{
public:
    virtual void           GetEnergy(EnergyBreakdown&,const DM_CD* cd ) const;
    virtual std::ostream&  Write    (std::ostream&) const;
    //! Whole-system RHF exchange via ERI4 bra-ket symmetry (doc/ERI4Rework.md §5.4), scaled -1/2 and
    //! cached per irrep; the 3-arg (energy) reuses the stashed whole basis; null basis -> per-irrep path.
    virtual const rsmat_t& GetMatrix(const obs_t*,const Spin&,const rChargeDensity*,const bs_t* wholeBasis) const;
    virtual const rsmat_t& GetMatrix(const obs_t*,const Spin&,const rChargeDensity*) const;
private:
    virtual rsmat_t CalcMatrix(const obs_t*,const Spin&,const rChargeDensity* cd) const;
    void EnsureWholeSystem(const rChargeDensity* cd) const;
    mutable size_t itsAllVersion=size_t(-1);
    mutable const bs_t* itsWholeBasis=nullptr;
    mutable std::map<std::string,rsmat_t> itsK;   //!< per-irrep exchange blocks, already scaled by -1/2
};

class VxcPol : public virtual Dynamic_HT, private Dynamic_HT_Imp_NoCache
{
public:
    virtual void           GetEnergy(EnergyBreakdown&,const DM_CD* cd ) const;
    virtual bool           IsPolarized() const {return true;}
    virtual std::ostream&  Write    (std::ostream&) const;
    //! Whole-system UHF exchange, PER SPIN (K^sigma from the sigma-density only), scaled -1 and cached
    //! per (spin,irrep); the 3-arg (energy) reuses the stash; null basis -> per-irrep path.
    virtual const rsmat_t& GetMatrix(const obs_t*,const Spin&,const rChargeDensity*,const bs_t* wholeBasis) const;
    virtual const rsmat_t& GetMatrix(const obs_t*,const Spin&,const rChargeDensity*) const;
private:
    virtual rsmat_t CalcMatrix(const obs_t*,const Spin&,const rChargeDensity* cd) const;
    void EnsureWholeSystem(const rChargeDensity* cd, const Spin& s) const;
    mutable size_t itsAllVersion=size_t(-1);
    mutable const bs_t* itsWholeBasis=nullptr;
    mutable std::map<Spin,std::map<std::string,rsmat_t>> itsK;   //!< per-spin, per-irrep K, scaled by -1
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
class FittedVee : public virtual Dynamic_HT, private Dynamic_HT_Imp
{
public:
    typedef std::shared_ptr<const BasisSet::FIT_CD_ABS> fbs_t;   //!< the charge-density (Coulomb-metric) fit face
    FittedVee(fbs_t& chargeDensityFitBasisSet, double numElectrons);
    ~FittedVee();   // anchored in the Imp TU (FittedCD complete there) so the unique_ptr can delete it
    virtual void          GetEnergy(EnergyBreakdown&,const DM_CD* cd) const;
    virtual std::ostream& Write    (std::ostream& os) const {return os;}
private:
    virtual rsmat_t CalcMatrix(const obs_t*,const Spin&,const rChargeDensity* cd) const;
    std::unique_ptr<ChargeDensity::FittedCD> itsFittedChargeDensity;   //!< owned (was a leaked raw ptr)
};


//###############################################################################
//
//  Linear least squares fit the unpolarized and polarized exchange-correlation potential.  The fit basis set is inserted by the
//  constructor and is not owned by FittedVxc.  The XC functional is owned by the inner LDAVxc.
//  Energy uses the exchange virial E_xc = 3/4 <rho|Vxc> -- correct for (pure) exchange; for correlation
//  use FittedVcorr below.
//
class FittedVxc : public virtual Dynamic_HT, private Dynamic_HT_Imp
{
public:
    typedef std::shared_ptr<const BasisSet::FIT_SF_ABS> fbs_t;   //!< the scalar-function (overlap-metric) fit face
    typedef std::shared_ptr<ExFunctional>     ex_t;

    FittedVxc(fbs_t& VxcFitBasisSet, ex_t&);
    ~FittedVxc();
    virtual void          GetEnergy       (EnergyBreakdown&,const DM_CD*) const;
    virtual void          UseChargeDensity(const rChargeDensity*);
    virtual std::ostream& Write           (std::ostream&) const;
private:
    virtual rsmat_t CalcMatrix(const obs_t*,const Spin&,const rChargeDensity*) const;

    std::unique_ptr<Fitting::FunctionFitter_Scalar<double>> itsFitter; //!< COMPOSED v_xc fit (was inherited)
    LDAVxc* itsLDAVxc;   //!< the v_xc=Vxc(rho) function to fit (concrete: it IS the ScalarFFClient)
};

class FittedVxcPol : public virtual Dynamic_HT, private Dynamic_HT_Imp_NoCache
{
public:
    typedef std::shared_ptr<const BasisSet::FIT_SF_ABS> fbs_t;   //!< the scalar-function (overlap-metric) fit face
    typedef std::shared_ptr<      ExFunctional>  ex_t;

    FittedVxcPol(fbs_t&, ex_t&);
   ~FittedVxcPol();
    // Required by HamiltonianTerm
    virtual void GetEnergy       (EnergyBreakdown&,const DM_CD* cd         ) const;
    virtual bool IsPolarized() const {return true;}

    virtual std::ostream&   Write(std::ostream&) const;
private:
    virtual rsmat_t CalcMatrix(const obs_t*,const Spin&,const rChargeDensity* cd) const;

    Dynamic_HT* itsUpVxc  ; //Spin up.
    Dynamic_HT* itsDownVxc; //Spin down.

};

//###############################################################################
//
//  This is a total energy term, not a matrix values Hamiltonian term.
//  A dedicated least-squares fit of an XC ENERGY DENSITY eps_xc(rho(r)) to the auxiliary fit basis,
//  exposed as a Dynamic_CC so the density contracts it:  E_xc = integral eps_xc rho = <rho|eps_xc_fit>.
//  Separate from a potential fit (different coefficients) but meant to SHARE its fit basis, so the
//  3-centre integrals are computed once.  Needed because eps_c != 3/4 v_c -- the exchange virial
//  (eps_x = 3/4 v_x) used by FittedVxc::GetEnergy is correct for exchange but wrong for correlation.
//  Composes a Fitting::FunctionFitter (from the Factory); Overlap is queried on it.
//
class FittedEpsXc : public virtual ChargeDensity::Dynamic_CC
{
public:
    typedef std::shared_ptr<const BasisSet::FIT_SF_ABS> fbs_t;   //!< the scalar-function (overlap-metric) fit face

    FittedEpsXc(fbs_t& fitBasisSet, const ExFunctional* ex);
    //! Re-fits eps_xc for this density and returns its matrix Sum_a c_a <Oi|f_a|Oj> for contraction.
    virtual const rsmat_t& GetMatrix(const obs_t*,const Spin&,const rChargeDensity* cd) const;
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
    virtual void GetEnergy(EnergyBreakdown&,const DM_CD* cd) const;
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
class FittedVcorrPol : public virtual Dynamic_HT, private Dynamic_HT_Imp_NoCache
{
public:
    typedef std::shared_ptr<const BasisSet::FIT_SF_ABS> fbs_t;   //!< the scalar-function (overlap-metric) fit face
    typedef std::shared_ptr<SpinCorrelation>            corr_t;  //!< the spin-native correlation functional

    FittedVcorrPol(fbs_t&, corr_t&);
   ~FittedVcorrPol();
    virtual void GetEnergy (EnergyBreakdown&, const DM_CD* cd) const;
    virtual bool IsPolarized() const {return true;}
    virtual std::ostream& Write(std::ostream&) const;
private:
    virtual rsmat_t CalcMatrix(const obs_t*, const Spin&, const rChargeDensity* cd) const;

    corr_t itsCorr;                                                      //!< the correlation functional (owned)
    std::unique_ptr<Fitting::FunctionFitter_Scalar<double>> itsVcFitter; //!< v_c^sigma potential fit
    std::unique_ptr<FittedEpsCPol>                          itsEpsC;     //!< E_c = integral eps_c rho (energy)
};


} //namespace