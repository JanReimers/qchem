// File: Hamiltonian/Internal/PWTerms.C  Plane-wave (dcmplx) Kohn-Sham Hamiltonian terms.
//
// These are the THIN terms that complete the dependency inversion: each derives from the dcmplx term
// base (cStatic_HT/cDynamic_HT in qcHamiltonian), holds the abstract orbital basis cobs_t, dynamic_casts
// it UP to the abstract BasisSet::Band_DFT_IBS<dcmplx> capability (in qcBasisSet), and asks that high-
// level question -- "the external matrix", "the Hartree matrix for this density".  The basis owns the
// integration; the term owns no G-vectors or mesh.  Energies delegate to the density's DM_Contract.
module;
#include <functional>
#include <iosfwd>
#include <memory>
export module qchem.Hamiltonian.Internal.PWTerms;
import qchem.Hamiltonian.Internal.Term;        // cStatic_HT / cDynamic_HT + their _Imp cache bases
import qchem.BasisSet.Band_FT_IBS;           // the reciprocal-space capability: Hartree/XC + external PP assembly
import qchem.BasisSet.Fit_IBS;               // cFIT_CD_ABS (the density-fit basis PW_Hartree is built with)
import qchem.Fitting.FunctionFitter;         // FunctionFitter_Density<dcmplx> (the fitter PW_Hartree holds, built once)
import qchem.Pseudopotential.Integrals_Pseudo;    // external-PP operator-assembly mixin + the local/separable models the term owns
import qchem.Hamiltonian.Internal.ExFunctional; // the LDA functional the XC term composes with the density
import qchem.Hamiltonian.Types;                 // cobs_t
import qchem.Structure;

export namespace qchem::Hamiltonian
{

//! External (pseudo)potential term for a plane-wave basis (static, density-independent).  THIS is the
//! pseudo-wall: the TERM owns the pseudopotential MODEL (an abstract local form factor + optional KB
//! nonlocal projector), and asks the basis to ASSEMBLE the matrix from it (MakeLocalPotential +
//! MakeSeparablePotential) -- physics lives Hamiltonian-side, integral assembly basis-side.  The models
//! are non-owning (the caller keeps them alive).  Pair with the kinetic, Hartree and XC terms for a full
//! Kohn-Sham Hamiltonian.
class PW_Pseudo
    : public virtual cStatic_HT
    , private        cStatic_HT_Imp
{
public:
    typedef std::shared_ptr<const Structure> st_t;
    PW_Pseudo(const st_t& st, const Pseudopotential::LocalPotential* loc,
                const Pseudopotential::SeparablePotential* nl=nullptr);
    virtual void          GetEnergy(EnergyBreakdown&, const cDM_CD*) const;
    virtual std::ostream& Write(std::ostream&) const;
private:
    virtual chmat_t CalculateMatrix(const cobs_t*, const Spin&) const;
    st_t theStructure;
    const Pseudopotential::LocalPotential*     itsLocal;       //!< local pseudopotential model (non-owning).
    const Pseudopotential::SeparablePotential* itsSep;         //!< KB nonlocal model (non-owning; may be null).
};

//! Non-relativistic kinetic ENERGY term T = 1/2 <p^2> for a plane-wave basis (diagonal in |k+G|^2).
//! Static (density-independent).  Uses the uncached MakeKinetic() block (the symmetric cache is bypassed
//! by the complex path).
class PW_Kinetic
    : public virtual cStatic_HT
    , private        cStatic_HT_Imp
{
public:
    virtual void          GetEnergy(EnergyBreakdown&, const cDM_CD*) const;
    virtual std::ostream& Write(std::ostream&) const;
private:
    virtual chmat_t CalculateMatrix(const cobs_t*, const Spin&) const;
};

//! Ion-ion (nuclear-nuclear) ENERGY term for a periodic crystal: the Ewald lattice sum of the ion cores
//! (charges = the cell atoms' itsZ = Zion).  The dcmplx sibling of Vnn -- it adds NO matrix contribution
//! (a constant), only the Enn energy, delegating to Ewald::NuclearRepulsion (which picks Ewald vs the
//! direct pair sum via Structure::isFinite()).
class PW_IonIon
    : public virtual cStatic_HT
    , private        cStatic_HT_Imp
{
public:
    typedef std::shared_ptr<const Structure> st_t;
    //! \a zionOf maps an atom's true species Z to its ION CORE charge (the PP's valence; identity for the
    //! all-electron baseline).  This is how the atom's itsZ stays the TRUE species while Ewald gets Zion.
    PW_IonIon(const st_t& st, std::function<double(int)> zionOf);
    virtual void          GetEnergy(EnergyBreakdown&, const cDM_CD*) const;
    virtual std::ostream& Write(std::ostream&) const;
private:
    virtual chmat_t CalculateMatrix(const cobs_t*, const Spin&) const;
    st_t theStructure;
    std::function<double(int)> itsZionOf;
};

//! Hartree (classical Coulomb) term for a plane-wave basis (density-dependent).  Asks the basis for the
//! Hartree matrix of the current density; the Poisson solve is the basis's business (G-space for PW).
class PW_Hartree
    : public virtual cDynamic_HT
    , private        cDynamic_HT_Imp
{
public:
    typedef std::shared_ptr<const BasisSet::cFIT_CD_ABS> fbs_t;
    //! Built with the density-fit basis obtained from the orbital basis's factory -- exactly as FittedVee
    //! takes its CD fit basis (BuildTerms creates it ONCE, never assuming orbital==fit).  The fitter is
    //! built here and reused every SCF cycle (DoFit the new density, then ask for the Hartree matrix).
    PW_Hartree(fbs_t chargeDensityFitBasisSet);
    ~PW_Hartree();
    virtual void          GetEnergy(EnergyBreakdown&, const cDM_CD*) const;
    virtual std::ostream& Write(std::ostream&) const;
private:
    virtual chmat_t CalcMatrix(const cobs_t*, const Spin&, const cChargeDensity*) const;
    std::unique_ptr<Fitting::FunctionFitter_Density<dcmplx>> itsFitter;   //!< ortho density fitter (built once)
};

//! Exchange-correlation term for a plane-wave basis, carrying ONE LDA functional (so a full LDA uses a
//! Dirac PW_XC + a VWN PW_XC, mirroring the molecular SlaterExchange+VWN split).  The matrix is the basis
//! integral of v_xc(rho(r)); the energy is integral eps_xc(rho) rho.  Both are real-space scalar fields
//! the term composes (functional o density) and hands to the basis -- the basis owns the integration.
class PW_XC
    : public virtual cDynamic_HT
    , private        cDynamic_HT_Imp
{
public:
    typedef std::shared_ptr<ExFunctional> xc_t;
    typedef std::shared_ptr<const BasisSet::cFIT_SF_ABS> fbs_t;
    //! Built with the Vxc fit basis obtained from the orbital basis's factory (BuildTerms creates it ONCE,
    //! never assuming orbital==fit) -- the overlap-metric sibling of PW_Hartree.
    PW_XC(const xc_t&, fbs_t vxcFitBasisSet);
    ~PW_XC();
    virtual void          GetEnergy(EnergyBreakdown&, const cDM_CD*) const;
    virtual std::ostream& Write(std::ostream&) const;
private:
    virtual chmat_t CalcMatrix(const cobs_t*, const Spin&, const cChargeDensity*) const;
    xc_t itsXc;
    std::unique_ptr<Fitting::FunctionFitter_Scalar<dcmplx>> itsScalarFitter;   //!< ortho scalar fitter (built once)
    //! The basis is captured from CalcMatrix so GetEnergy (which has no basis parameter) can ask it for
    //! the energy integral integral eps_xc rho with the current density.  Same basis every iteration.
    mutable const BasisSet::Band_FT_IBS* itsBasis=nullptr;
};

} //namespace
