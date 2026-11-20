// Hamiltonian/Internal/Terms.C  Declare and export all Hamiltonian term types.
module;
#include <iosfwd>
#include <memory>
#include <type_traits>
#include <map>
#include <string>
#include <vector>
#include <functional>
export module qchem.Hamiltonian.Internal.Terms;
import qchem.Hamiltonian.Internal.Term;
import qchem.Hamiltonian.Internal.GeometryGradient;   // Ven's force face (GeometryOptimizationPlan §13.2; lane B)
import qchem.Hamiltonian.Internal.ExFunctional;
import qchem.Structure;
import qchem.Fitting.FunctionFitter;          // Fitting::FunctionFitter (composed; clients never see the impl)
import qchem.Fitting.CoulombOperator;         // Fitting::FirstOrderDensityFit (what Vee_Fit holds)
import qchem.ChargeDensity;
import qchem.ChargeDensity.TransitionDensity;   // δD: the response face (LinearResponsePlan C1/H2)
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
          , public virtual tGeometryGradient_HT<double>, public virtual tStatic_PC<double>   // its force (§13.2; lane B): Imp/Ven_Forces.C
{
public:
    typedef std::shared_ptr<const Structure> st_t;
    Ven(const st_t& st);
    virtual void          GetEnergy(EnergyBreakdown&,const rDM_CD* cd) const override;
    virtual std::ostream& Write    (std::ostream&) const override;
    // The force (GeometryOptimizationPlan §13.2; lane B): Imp/Ven_Forces.C.
    virtual std::string MissingGradient(const tbs_t<double>* bs, const Structure* st) const override;
    virtual void        AddForces(rvec3vec_t& F, const rDM_CD* D, const Structure* st) const override;
    virtual rvec3vec_t  Pullback(const robs_t* bs, const Spin& s, const hmat_t<double>& D, const Structure* st) const override;
private:
    virtual rsmat_t MakeMatrix(const robs_t*,const Spin&) const override;
    st_t theStructure;
};

//###############################################################################
//
//  Local pseudopotential electron-ion term: the pseudized replacement for Ven.  Instead of the analytic
//  -Z/r nuclear attraction: a basis that answers Orbital_PP_IBS assembles V_loc in closed form; otherwise it
//  quadratures the smooth real-space V_loc(r) on the molecular/atomic mesh,
//  <chi_i|V_loc|chi_j> = Sum_g w_g chi_i(r_g) chi_j(r_g) V_loc(r_g) (= the XC-path MatrixOverlap shape).
//  STATIC (density-independent), so it is built once.  V_loc is the species radial field's real-space view (SpeciesRadialField::ValueR).
//
class PP_Local : public virtual rStatic_HT, private rStatic_HT_Imp
              , public virtual tGeometryGradient_HT<double>, public virtual tStatic_PC<double>   // its force (§13.2; lane B): Imp/PP_Forces.C
{
public:
    // The force (GeometryOptimizationPlan §13.2; lane B): Imp/PP_Forces.C.
    virtual std::string MissingGradient(const tbs_t<double>* bs, const Structure* st) const override;
    virtual void        AddForces(rvec3vec_t& F, const rDM_CD* D, const Structure* st) const override;
    virtual rvec3vec_t  Pullback(const robs_t* bs, const Spin& s, const hmat_t<double>& D, const Structure* st) const override;
    //! The virial theorem needs a Coulombic (degree -1 homogeneous) potential; a pseudopotential is not
    //! (erf-screened local part + KB projectors), so the SCF drops both the virial gate and column (V1.27).
    virtual bool IsVirialValid() const override {return false;}
    typedef std::shared_ptr<const Structure> st_t;
    typedef std::shared_ptr<const BasisSet::SpeciesRadialField> vloc_t;
    typedef std::shared_ptr<const qcMesh::Mesh> mesh_t;
    //! \a mesh: the run's PP quadrature, built ONCE by the Hamiltonian in its own build phase and shared with PP_NonLocal (it is
    //! geometry-only, so one mesh serves every basis block) -- never built lazily at the first matrix request.  Null when every
    //! basis block assembles the field in closed form (\c Orbital_PP_IBS).
    PP_Local(const st_t& st, vloc_t vloc, mesh_t mesh);
    virtual void          GetEnergy(EnergyBreakdown&,const rDM_CD* cd) const override;   // Een (PP local) = DM_Contract
    virtual std::ostream& Write    (std::ostream&) const override;
private:
    virtual rsmat_t MakeMatrix(const robs_t*,const Spin&) const override;
    st_t             theStructure;
    vloc_t           itsVloc;
    mesh_t           itsMesh;
};

//###############################################################################
//
//  Separable (Kleinman-Bylander) NON-LOCAL pseudopotential term.  V_NL = Sum_{a,p} D_p Sum_m |beta_p Y_lm>
//  <beta_p Y_lm| -- real symmetric, STATIC.  This is the repulsive (for the occupied valence l-channels)
//  piece that lifts the over-bound local-only spectrum back to the all-electron valence eigenvalues.
//
//  THREE assemblies, chosen by capability cross-casts on the block:
//   * CLOSED FORM (a basis that answers Orbital_PP_IBS, the molecular Gaussians): exact, no mesh.
//   * EXPLICIT-angular (molecular/Cartesian) -- the 3-D mesh route: per atom, per projector, per m=-l..l a
//     rank-1 D|b><b| with b_i = <chi_i|beta_p(|r-R|) Y_lm(rhat)> (mesh quadrature).
//   * IMPLICIT-angular (ATOMIC, radial: BasisSet::ImplicitAngular_IBS) -- the per-l RADIAL route, because a
//     radial block's stored chi_i OMIT the irrep's Y_lm and the 3-D route would silently give nonsense
//     (an l=0 projector leaking into every l block; every l>=1 projector integrating to zero -- the
//     occupied-d defect found via MnO, doc/SymmetryUpgradePlan.md sec 7 step 7).  See MakeMatrixRadial.
//  Verified against the reciprocal (2l+1)P_l form; per-l gate A_PP.PerLKleinmanBylanderOracle (s,p,d,f).
//
class PP_NonLocal : public virtual rStatic_HT, private rStatic_HT_Imp
                 , public virtual tGeometryGradient_HT<double>, public virtual tStatic_PC<double>   // its force (§13.2; lane B): Imp/PP_Forces.C
{
public:
    // The force (GeometryOptimizationPlan §13.2; lane B): Imp/PP_Forces.C.
    virtual std::string MissingGradient(const tbs_t<double>* bs, const Structure* st) const override;
    virtual void        AddForces(rvec3vec_t& F, const rDM_CD* D, const Structure* st) const override;
    virtual rvec3vec_t  Pullback(const robs_t* bs, const Spin& s, const hmat_t<double>& D, const Structure* st) const override;
    //! The virial theorem needs a Coulombic (degree -1 homogeneous) potential; a pseudopotential is not
    //! (erf-screened local part + KB projectors), so the SCF drops both the virial gate and column (V1.27).
    virtual bool IsVirialValid() const override {return false;}
    typedef std::shared_ptr<const Structure> st_t;
    typedef std::shared_ptr<const BasisSet::SpeciesProjectorSet_R> sep_t;
    typedef std::shared_ptr<const qcMesh::Mesh> mesh_t;
    //! \a mesh: the run's shared PP quadrature (see PP_Local) for the explicit-angular route, null when no block needs it; \a mp still states the radial
    //! route's radial rule (an atomic block never touches the 3-D mesh).
    PP_NonLocal(const st_t& st, sep_t sep, const qcMesh::MeshParams& mp, mesh_t mesh);
    virtual void          GetEnergy(EnergyBreakdown&,const rDM_CD* cd) const override;   // Een (PP nonlocal) = DM_Contract
    virtual std::ostream& Write    (std::ostream&) const override;
private:
    virtual rsmat_t MakeMatrix(const robs_t*,const Spin&) const override;
    //! The ATOMIC (implicit-angular / radial) assembly -- see the body: a radial block's stored functions
    //! omit the irrep's Y_lm, so the 3-D mesh route cannot form <chi|beta Y_lm> and must not be used.
    rsmat_t MakeMatrixRadial(const BasisSet::ImplicitAngular_IBS&, size_t n) const;
    st_t             theStructure;
    sep_t            itsSep;
    qcMesh::MeshParams itsMeshParams;   //!< the radial route's radial rule
    mesh_t           itsMesh;           //!< the explicit-angular route's shared PP quadrature
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
class Dynamic_HF_HT_Imp
    : public virtual rDynamic_HF_HT
    , public virtual tResponse_HT<double>   //!< J and K are LINEAR in D: their response is their own sweep of δD
{
public:
    //! Fock build: assemble the whole-system blocks ONCE per (density, spin) from the composite \a wholeBasis
    //! using ERI4 bra-ket symmetry (canonical pairs -> ScatterBoth), cache the per-irrep blocks, return this
    //! irrep's block.  \a wholeBasis is required (HF is whole-system); a null basis throws.
    virtual const rsmat_t& GetMatrix(const robs_t*,const Spin&,DensityID,const rbs_t* wholeBasis) const override;
    //! The ground-state refresh: latch \a wholeBasis (null = keep the latched one) and build the blocks of every spin irrep
    //! \a cd resolves (\c SpinIrrepsOf, folded through \c CacheSpin).
    virtual void RefreshForDensity(const rbs_t* wholeBasis, const rChargeDensity* cd) const override;
    //! \name THE RESPONSE (doc/LinearResponsePlan.md H2): J[δD] / K[δD] through the SAME contraction body,
    //! into a SEPARATE cache slot keyed on δ's serial -- a response build never overwrites the ground-state
    //! blocks.  \a D0 is not needed: the operators are linear.  q = 0 only (bra == ket).
    //!@{
    virtual void   RefreshForDensity(const rbs_t* wholeBasis, const rChargeDensity* D0,
                                     const ChargeDensity::TransitionDensity<double>& delta) const override;
    virtual rmat_t GetMatrix(const robs_t* bra, const robs_t* ket, const Spin& s,
                             const ChargeDensity::TransitionDensity<double>& delta) const override;
    //!@}
protected:
    //! The one operation that distinguishes Coulomb from exchange: scatter the \a sweep operand across the
    //! canonical irrep pairs into the zeroed per-irrep blocks \a X (one per irrep, in the operand's walk order).
    //! It takes the narrow HF face -- the ONLY thing the scatter consumes -- so the SAME call serves a ground-
    //! state density and a transition density (doc/LinearResponsePlan.md §3c: the J/K chain is not duplicated).
    virtual void   AccumulateAll(std::vector<rsmat_t>& X,const ChargeDensity::tHF_System_CD<double>& sweep) const=0;
    //! \name THE SPIN AXIS (V1.37 step 3) -- two questions, answered by the concrete term, that used to
    //! be answered by WHICH TYPE it was (Vee / Vxc(-1/2) / two Vxc(-1) inside a VxcPol).
    //!@{
    //! Which spin keys my blocks, and so which channel of the density they are built FROM: \c Spin::None when
    //! the operator ignores the channel (Coulomb sees the total), \a s itself when it is per channel (exchange
    //! is same-spin).
    virtual Spin           CacheSpin (const Spin& s) const=0;
    //! Fock coefficient applied to every block after the scatter (1 for Coulomb; the K coefficient for Vxc,
    //! which depends on whether the block is a folded doublet or one channel).
    virtual double         Scale(const Spin& s) const=0;
    //!@}
    //! The density my spin-\a s blocks are built from: the \c CacheSpin(s) channel of \a cd (the whole of \a cd
    //! for \c Spin::None).  THROWS if \a cd does not resolve that channel or carries no density matrix.
    const rDM_CD*  DensityFor(const rChargeDensity* cd, const Spin& s) const;
    //! The blocks the refresh built for \a cd and spin \a s (a CacheSpin): a READ, for \c GetMatrix and \c GetEnergy.  THROWS if
    //! the refresh never saw \a cd (a caller that skipped \c RefreshForDensity).
    const std::map<Irrep,rsmat_t>& RefreshedBlocks(DensityID id, const Spin& s) const;

    mutable const rbs_t* itsWholeBasis=nullptr;    //!< whole basis (latched by the first Fock refresh; stable across the run)
    //! The J (Coulomb) or K (exchange) per-irrep blocks of ONE spin: accumulated (over irreps) and contracted
    //! (over that spin's D) for the density ID'd by \c version.  Keyed by the block's spatial Irrep
    //! (\c GetIrrep(Spin::None)), already scaled.
    struct Blocks
    {
        size_t version=size_t(-1);   //!< the operand's serial (a ground-state density, or a transition density for the response)
        DensityID id=DensityID(static_cast<const rChargeDensity*>(nullptr));   //!< a ground-state density's identity (GetMatrix's check)
        std::map<Irrep,rsmat_t> jk;
    };
    mutable std::map<Spin,Blocks> itsJKs;          //!< one Blocks per CacheSpin this run asks for
private:
    //! Latch the composite basis on first use; THROW on a changed one, or on a null one before any latch.
    void LatchWholeBasis(const rbs_t* wholeBasis) const;
    //! Contract \a cd into the whole-system blocks for spin \a s (a CacheSpin) if stale for this density (the refresh's body).
    const std::map<Irrep,rsmat_t>& ContractAll(const rChargeDensity* cd, const Spin& s) const;
    //! The response blocks of CacheSpin \a s for \a delta (its \a s channel through the HF sweep face).
    const std::map<Irrep,rsmat_t>& ContractResponse(const ChargeDensity::TransitionDensity<double>& delta,
                                                          const Spin& s) const;
    mutable std::map<Spin,Blocks> itsResponseJKs;   //!< the response twin of itsJKs, keyed on δ's serial
    //! THE ONE CONTRACTION BODY: scatter \a sweep (identified by \a version) into \a slot for spin \a s, unless
    //! \a slot already holds that version.  Every entry -- the Fock build, the energy, a response -- comes here.
    const std::map<Irrep,rsmat_t>& Contract(Blocks& slot, const ChargeDensity::tHF_System_CD<double>& sweep,
                                                  size_t version, const Spin& s) const;
};

class Vee_HF : public Dynamic_HF_HT_Imp
         , public virtual tGeometryGradient_HT<double>   // its force (§13.2; lane B): Imp/HF_Forces.C
{
public:
    // The force (GeometryOptimizationPlan §13.2; lane B): Imp/HF_Forces.C.
    virtual std::string MissingGradient(const tbs_t<double>* bs, const Structure* st) const override;
    virtual void        AddForces(rvec3vec_t& F, const rDM_CD* D, const Structure* st) const override;
    virtual void          GetEnergy(EnergyBreakdown&,const rDM_CD* cd ) const override;
    virtual std::ostream& Write    (std::ostream&) const override;
protected:
    virtual void          AccumulateAll(std::vector<rsmat_t>& X,const ChargeDensity::tHF_System_CD<double>& sweep) const override;
    // Coulomb sees the TOTAL density and ignores the channel: one set of blocks serves every spin.
    virtual Spin          CacheSpin (const Spin&) const override {return Spin::None;}
    virtual double        Scale(const Spin&) const override {return 1.0;}
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
class Vxc_HF : public Dynamic_HF_HT_Imp
         , public virtual tGeometryGradient_HT<double>   // its force (§13.2; lane B): Imp/HF_Forces.C
{
public:
    // The force (GeometryOptimizationPlan §13.2; lane B): Imp/HF_Forces.C.
    virtual std::string MissingGradient(const tbs_t<double>* bs, const Structure* st) const override;
    virtual void        AddForces(rvec3vec_t& F, const rDM_CD* D, const Structure* st) const override;
    virtual void           GetEnergy(EnergyBreakdown&,const rDM_CD* cd ) const override;
    virtual std::ostream&  Write    (std::ostream&) const override;
protected:
    virtual void          AccumulateAll(std::vector<rsmat_t>& X,const ChargeDensity::tHF_System_CD<double>& sweep) const override;
    virtual Spin          CacheSpin (const Spin& s) const override {return s;}          // same-spin: per channel
    virtual double        Scale(const Spin& s) const override {return s==Spin::None ? -0.5 : -1.0;}
};

//###############################################################################
//
//  Implementation of the Coulomb potential
//
//            /
// Vee(r_1) = | Ro_fit(r_2)/r_12 d^3 r_2
//           /
//
// Where Ro is actually a fitted charge density.  Two terms, because two formulas (V1.1b, ruled 2026-11-06; shapes in
// doc/Records/HartreeFitConsistency.md).  Neither serves a solid: neither constructor takes a complex fit, and the periodic route is
// Vee_Fit<dcmplx>.
//
// DUNLAP, C^T G^-1 C.  Ro_fit is expanded in an auxilliary basis in the Coulomb metric, charge constrained, so the matrix elements
// are three center integrals and the four center ones of Hartree-Fock never appear.  Second order in the fit error.  Refused on
// solids by construction: a periodic Gaussian fit has no (ij|f) to call, and on plane waves this and Vee_Fit are one formula.
//
class Vee_FitDunlap : public virtual rDynamic_HT, private rDynamic_HT_Imp
                    , public virtual tGeometryGradient_HT<double>, public virtual tStatic_PC<double>   // its force (§13.2; lane B): Imp/Vee_FitDunlap_Forces.C
{
public:
    // The force (GeometryOptimizationPlan §13.2; lane B): Imp/Vee_DFT_Dunlap_Forces.C.
    virtual std::string MissingGradient(const tbs_t<double>* bs, const Structure* st) const override;
    virtual void        AddForces(rvec3vec_t& F, const rDM_CD* D, const Structure* st) const override;
    virtual rvec3vec_t  Pullback(const robs_t* bs, const Spin& s, const hmat_t<double>& D, const Structure* st) const override;
    typedef std::shared_ptr<const BasisSet::rFIT_CD_ABS> fbs_t;   //!< the charge-density (Coulomb-metric) fit face
    typedef std::unique_ptr<Fitting::DunlapDensityFit<double>> FitPtr;
    //! The fitter is the factory's choice (\c Fitting::MakeDunlapDensityFitter); this term holds only its face (ruling (g)).
    //! \a numElectrons is the FIXED N the fit is constrained to (user ruling 2026-11-08, via lane B: not the density's Tr(DS),
    //! which moves with the geometry at fixed D and would put a λ Tr(D ∂S/∂R) term into the force).
    Vee_FitDunlap(FitPtr fit, double numElectrons);
    ~Vee_FitDunlap();
    //! \copydoc tDynamic_HT::RefreshForDensity
    //! The fitted charge density is k-INDEPENDENT and refit once per density serial, so it is exactly what
    //! this phase is for.  It used to be fit LAZILY inside \c MakeMatrix, i.e. from inside the block loop.
    virtual void RefreshForDensity(const rChargeDensity* cd) const override;
    //! \f$2E_{fit}-E_{fitfit}\f$, which needs \a cd's fit to be the Coulomb one.  A seed (\c NumericCD) is overlap-fitted, and the formula
    //! would be wrong for it: taking a density MATRIX here is the one thing keeping a seed out, so this argument type is load-bearing.
    virtual void          GetEnergy(EnergyBreakdown&,const rDM_CD* cd) const override;
    virtual std::ostream& Write    (std::ostream& os) const override {return os;}
    //! Hartree depends on the TOTAL density: one block for every spin channel.
    virtual Spin CacheSpin(const Spin&) const override {return Spin::None;}
    //! The ENERGY block GetEnergy contracts: BUILT from the fit the energy pass refreshed, with no memo -- the energy refresh has
    //! no basis, so it builds no blocks (V1.41 stage 3; the same shape as FittedVxc's eps block).
    virtual const rsmat_t& GetEMatrix(const robs_t* bs,const Spin& s,DensityID id) const override
    {itsEMat=MakeMatrix(bs,s,id); return itsEMat;}
private:
    mutable rsmat_t itsEMat;   //!< GetEMatrix's returned block (one serial contraction reads it at once)
    virtual rsmat_t MakeMatrix(const robs_t*,const Spin&,DensityID) const override;
    //! The density projects ITSELF onto the fitter, and the fit is constrained to the fixed electron count N.
    void Fit(const rChargeDensity& cd) const {itsFit->DoFit(cd.Project(*itsFit), itsNumElectrons);}
    FitPtr itsFit;   //!< owned
    double itsNumElectrons;   //!< the fixed N of the charge constraint
};

//
// NORMAL, O^T K O.  q = O D, v = K q, H = O^T v: rho on a mesh, K a CoulombOperator, whose energy this term asks for.  First order in
// the mesh.  Unconstrained, so it reports the mesh's charge error (\c HartreeChargeError) and the first-order energy effect of it
// (\c HartreeChargeEnergy, \f$2E_H\,\delta q/N\f$; ruled 2026-11-06, reported not rescaled).
//
//! \brief THE RESPONSE HALF of the periodic fitted Hartree term (the user's \c Vee_Response, 2026-11-08): the transition density's
//! own G+q Poisson, \f$\delta V_H(G+q)=4\pi\,\delta\tilde\rho(G+q)/|G+q|^2\f$, then each (k+q, k) pair's exact adjoint
//! (doc/LinearResponsePlan.md R2, R3).  Linear, so no D0.  A Bloch-q field (\c ΔGq_Map), so it stays its own machinery until the
//! response route is redone; \c Vee_Fit<dcmplx> carries it as a mixin.
class Vee_Response : public virtual tResponse_HT<dcmplx>
{
public:
    explicit Vee_Response(std::shared_ptr<const BasisSet::cFIT_CD_ABS> fb) : itsFitBasis(std::move(fb)) {}
    virtual void RefreshForDensity(const cbs_t* wholeBasis, const cChargeDensity* D0, const TransitionDensity<dcmplx>& delta) const override;
    virtual mat_t<dcmplx> GetMatrix(const cobs_t* bra, const cobs_t* ket, const Spin& s, const TransitionDensity<dcmplx>& delta) const override;
private:
    std::shared_ptr<const BasisSet::cFIT_CD_ABS> itsFitBasis;   //!< the CD fit basis the transition density projects on
    mutable size_t  itsDeltaVersion=size_t(-1);                 //!< the transition density's serial \c itsDeltaField holds
    mutable ΔGq_Map itsDeltaField;                              //!< \f$\delta V_H(G+q)\f$ of that transition density
};

//! The per-run-scalar corners of \c Vee_Fit: nothing on the molecular instantiation; on the periodic one, the real TRIM block's
//! cache (its matrix is the same \c MakeBlock) and the response half.
template <class T, class D> class VeeFitCorners
{
public:
    explicit VeeFitCorners(std::shared_ptr<const BasisSet::FIT_CD_ABS<T>>) {}
};
template <class D> class VeeFitCorners<dcmplx,D>
    : public Dynamic_HT_RealBlock_Imp
    , public Vee_Response
{
public:
    explicit VeeFitCorners(std::shared_ptr<const BasisSet::cFIT_CD_ABS> fb) : Vee_Response(std::move(fb)) {}
    //! The real block's ENERGY block: built, as the scalar one is (the energy refresh builds no blocks).
    virtual const rsmat_t& GetEMatrix(const robs_t* bs,const Spin& s,DensityID id) const override
    {itsEMatR=MakeMatrix(bs,s,id); return itsEMatR;}
protected:
    virtual rsmat_t MakeMatrix(const robs_t* bs,const Spin& s,DensityID id) const override
    {return static_cast<const D&>(*this).template MakeBlock<double>(bs,s,id);}
private:
    mutable rsmat_t itsEMatR;
};

//
// THE FIRST-ORDER FITTED HARTREE TERM, GENERIC OVER THE FIT BASIS (V-ONE-OF-N 3b; user rulings 2026-11-08 (c), (e), (f), (g),
// doc/Records/VOneOfN_D2Design.md §1, §1a).  It was Vee_DFT (molecular δ mesh) and Vee_Hartree (plane wave / GPW), the same
// formula written twice.  ONE code path, no basis-kind branch:
//   - the density projects itself onto the fitter (\c tChargeDensity::Project), and the fitter keeps its coefficients;
//   - the matrix is the PHYSICS round: the fitter contracts them with its basis's COULOMB 3-centre integrals,
//     \f$\sum_f c_f\langle ab|1/r|f\rangle\f$ -- a δ mesh through its Coulomb operator, a plane wave as \f$k(G)\langle ab|G\rangle\f$;
//   - the energy is the fitter's first-order \f$E_H\f$; \f$\mathrm{Tr}(DV)\f$ is computed only for a fitter that reads it.
// Unconstrained, so it reports the fit's charge error (\c HartreeChargeError) and its first-order energy effect
// (\c HartreeChargeEnergy, \f$2E_H\,\delta q/N\f$; ruled 2026-11-06, reported not rescaled), when \a numElectrons is known (> 0).
//
template <class T> class Vee_Fit
    : public virtual tDynamic_HT<T>
    , private        tDynamic_HT_Imp<T>
    , public         VeeFitCorners<T, Vee_Fit<T>>
    , public virtual tGeometryGradient_HT<T>            // its force (§13.2; lane B 2026-11-09): Imp/Vee_Fit_Forces.C
    , public virtual tStatic_PC<T>                      // the per-block pullback client DM_Pullback walks
    , public         RealPullbackBase<T,Vee_Fit<T>>     // ... and on <dcmplx> the real TRIM block's
    , public virtual tStrainGradient_HT<T>              // its virial (CellHalf_Design §5; lane H 2026-11-10): Imp/Vee_Fit_Virials.C
    , public virtual tStatic_VC<T>                      // the per-block virial client DM_Virial walks
    , public         RealVirialBase<T,Vee_Fit<T>>       // ... and on <dcmplx> the real TRIM block's
{
public:
    // The virial (CellHalf_Design §5): the fitter's FIT side (MakeHartreeVirial) plus every block's ORBITAL side (MakeRepulsionVirial).
    virtual std::string MissingVirial(const tbs_t<T>* bs, const Structure* st) const override;
    virtual void        AddVirial(rmat3d_t& V, const tDM_CD<T>* D, const Structure* st) const override;
    virtual rmat3d_t    Virial(const tobs_t<T>* bs, const Spin& s, const hmat_t<T>& D, const Structure* st) const override
    {return VirialBlock<T>(bs, s, D, st);}
    //! ONE virial body for both block scalars, as PullbackBlock.  Public for the mixin.
    template <class U> rmat3d_t VirialBlock(const tobs_t<U>* bs, const Spin&, const hmat_t<U>& D, const Structure* st) const;
    // The force (GeometryOptimizationPlan §13.2; lane B): Imp/Vee_Fit_Forces.C.
    virtual std::string MissingGradient(const tbs_t<T>* bs, const Structure* st) const override;
    virtual void        AddForces(rvec3vec_t& F, const tDM_CD<T>* D, const Structure* st) const override;
    virtual rvec3vec_t  Pullback(const tobs_t<T>* bs, const Spin& s, const hmat_t<T>& D, const Structure* st) const override
    {return PullbackBlock<T>(bs, s, D, st);}
    //! ONE pullback body for both block scalars (the native block, and on <dcmplx> the real TRIM block): the fitter's.  Public for the mixin.
    template <class U> rvec3vec_t PullbackBlock(const tobs_t<U>* bs, const Spin&, const hmat_t<U>& D, const Structure* st) const;
    using FitPtr = std::unique_ptr<Fitting::FirstOrderDensityFit<T>>;
    //! \a responseFit: the CD fit basis the periodic response half projects δρ on (unused, and may be null, on a molecule).
    Vee_Fit(FitPtr fit, double numElectrons, std::shared_ptr<const BasisSet::FIT_CD_ABS<T>> responseFit=nullptr);
    ~Vee_Fit();
    using tDynamic_HT<T>::RefreshForDensity;
    //! \copydoc tDynamic_HT::RefreshForDensity
    virtual void RefreshForDensity(const tChargeDensity<T>* cd) const override;
    //! \copydoc HT_BlockOwner::BuildBlocks -- on the periodic instantiation I own two caches (Bloch, real TRIM) and prepare both.
    virtual void BuildBlocks(const tbs_t<T>* bs, const tChargeDensity<T>* cd) const override
    {
        tDynamic_HT_Imp<T>::BuildBlocks(bs, cd);
        if constexpr (std::is_same_v<T,dcmplx>) Dynamic_HT_RealBlock_Imp::BuildBlocks(bs, cd);
    }
    virtual void          GetEnergy(EnergyBreakdown&,const tDM_CD<T>* cd) const override;
    virtual std::ostream& Write    (std::ostream& os) const override;
    //! Hartree depends on the TOTAL density: one block for every spin channel.
    virtual Spin CacheSpin(const Spin&) const override {return Spin::None;}
    //! The ENERGY block (read only by a fitter whose energy needs Tr(DV), and by the eh_trace A/B): BUILT from the fit the energy pass
    //! refreshed, with no memo -- the energy refresh has no basis, so it builds no blocks (V1.41 stage 3).
    virtual const hmat_t<T>& GetEMatrix(const tobs_t<T>* bs,const Spin& s,DensityID id) const override
    {itsEMat=MakeMatrix(bs,s,id); return itsEMat;}
    //! One block's matrix, for either block scalar \a U on this run's fit axis: the fitter's physics round.  Public for the corner.
    template <class U> hmat_t<U> MakeBlock(const tobs_t<U>* bs, const Spin&, DensityID id) const;
private:
    mutable hmat_t<T> itsEMat;   //!< GetEMatrix's returned block (one serial contraction reads it at once)
    virtual hmat_t<T> MakeMatrix(const tobs_t<T>* bs,const Spin& s,DensityID id) const override {return MakeBlock<T>(bs,s,id);}
    void Fit(const tChargeDensity<T>& cd) const;
    FitPtr itsFit;
    double itsNumElectrons;
};


// FittedVxc (the least-squares v and eps fits of the A1 / multipole cells) is RETIRED (lane S, 2026-11-15): those cells run the one XC term,
// tVxc_DFT, over the fitted-potential sampler (ChargeDensity::MakeFittedPotentialSampler / MakeMultipoleFittedPotentialSampler).

} //namespace