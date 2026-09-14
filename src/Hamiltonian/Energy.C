// File: Hamiltonian/Energy.C  The energy (and charge) accounting of one density: keyed contributions with ROLES.
//
// V1.12 (user review 2026-09-13).  This used to be a struct of 13 public doubles, of which eight summed to the
// total, five were diagnostics that must NOT be summed (EenNL is a subset of Een; the Dunlap fit pieces are
// already combined into Eee) and one was not an energy at all.  Every new term family edited the struct, the
// totals, op+= and Display.  Now a term INSERTS its contribution under a unique name with a ROLE, and the
// totals are role sums -- a relativistic Hamiltonian adds "RestMass", a DFT one "Exc" (or "Eex"+"Ecorr"),
// +U adds "E_U", and nothing here changes.
module;
#include <map>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>
export module qchem.Energy;
import qchem.Types;   // rvec_t -- the per-site moments

export namespace qchem
{

//! \brief WHAT KIND of energy a contribution is -- the axis the totals are summed on.  A ROLE, not a name
//! prefix: a naming convention is one nobody compiles.
enum class EnergyRole
{
    Kinetic,    //!< \f$\langle T\rangle\f$ (NR \f$\tfrac12\langle p^2\rangle\f$, or the Dirac kinetic energy)
    Potential,  //!< an electronic potential energy: electron-ion, Hartree, exchange, correlation, +U, ...
    Constant,   //!< a density-INDEPENDENT constant of the structure: the ion-ion Madelung/Coulomb \f$E_{nn}\f$,
                //!< the periodic \f$G=0\f$ alignment \f$E_{\alpha Z}\f$ (QE/CP2K's "alpha Z" -- the finite remainder
                //!< of the divergent electron/ion \f$G=0\f$ terms once the neutralising background is imposed)
    RestMass,   //!< \f$\langle mc^2\rangle\f$ of a Dirac run
    Entropy     //!< the Mermin \f$-TS\le0\f$ of Fermi smearing: the ONLY footprint of smearing on the energy (the
                //!< entropy never touches \f$H\f$; \f$f\f$ enters only \f$D\f$).  Stamped by the iterator at fill time.
};

//! \brief One contribution.  \c TrDV is the term's EXPECTATION in the Fock operator, \f$\mathrm{Tr}(D\,V_{\rm term})\f$
//! -- the second number the BAND-ENERGY form of the total needs (\c GetBandEnergy).  Optional: a term supplies it
//! where it is free (a linear term: \f$=E\f$; a Coulomb/exchange quadratic: \f$2E\f$; a constant: 0) and leaves it
//! absent where it is not (the fitted Hartree and the XC terms under density MIXING, where \f$\mathrm{Tr}(D_{out}
//! V[\rho_{in}])\f$ is the Harris–Foulkes subtlety and is not to be guessed).  Absent is honest; the band form
//! then throws naming the term, rather than silently evaluating a different functional.
struct EnergyTerm
{
    double                E    = 0.0;
    std::optional<double> TrDV;
    EnergyRole            role = EnergyRole::Potential;
};

//! \brief The CHARGE accounting of the same density -- the structure the user asked for
//! (\f$\{N,\rho_\uparrow,\rho_\downarrow,\text{lost},\text{atoms}\{\rho_{i\uparrow},\rho_{i\downarrow},m_i\}\}\f$), grown as far
//! as the terms can fill it.  \c lost is the signed charge the collocation grid could not represent,
//! \f$\int\tilde\rho\,d^3r-\mathrm{Tr}(DS)\f$ (== CP2K's "Electronic density on regular grids" error; 0 on a
//! gridless path).  "Grid" is the MECHANISM, so it is not in the name.
//!
//! \c siteMoments is THE OBSERVABLE OWNER for the integrated per-site spin moments (R1.0h, 2026-09-14):
//! \f$\mu_A=\int w_A(\rho_\uparrow-\rho_\downarrow)\,d^3r\f$ in electrons, one entry per site block of the
//! XC quadrature's atom-centred partition, EMPTY when the run has none (a uniform raster, or an unpolarized
//! run -- both "not measurable", never "zero").  Filled by the spin-native XC term in its ENERGY pass, where
//! both channel rasters are already in hand -- so the number rides the same \c EnergyBreakdown the SCF
//! trace and the observer already receive, contemporaneously with the iteration that produced it.  (It
//! used to be PULLED through three faces -- Hamiltonian -> term -> sampler -- and reported from inside the
//! sampler's cache-advance branch; a value channel replaced the pull, the trace replaced the reach-in.)
//! Named partition (Becke), because until Bader's zero-flux basins land the number is partition-dependent.
struct ChargeBreakdown
{
    double lost = 0.0;
    rvec_t siteMoments;
};

//! \brief The energy breakdown of one density: keyed, role-tagged, insertion-ORDERED contributions (Display
//! reads in term order), plus a second map of DIAGNOSTICS that are reported but never summed.
class EnergyBreakdown
{
public:
    //! A term's contribution.  Merges by name (\f$+=\f$): the two spin channels of a polarized term, or the
    //! per-irrep pieces of one, land in ONE entry.  The role must agree with an existing entry's.
    void Add(const std::string& name, double E, EnergyRole role, std::optional<double> TrDV = std::nullopt);
    //! A diagnostic beside the contributions: reported, NEVER summed (a sub-split such as the nonlocal part of
    //! \f$E_{en}\f$, or the Dunlap fit pieces whose combination is already in the Hartree entry).  Merges (+=).
    void AddDiagnostic(const std::string& name, double v);

    double operator[](std::string_view name) const;   //!< a contribution's E; 0 if absent
    double Diagnostic(std::string_view name) const;   //!< a diagnostic;       0 if absent
    bool   Has(std::string_view name) const;

    //! Total energy = the sum of EVERY contribution.  With Fermi smearing this is the Mermin FREE ENERGY
    //! \f$A=E-TS\f$ -- the quantity the finite-T SCF makes stationary -- because the Entropy entry is in the sum;
    //! with no smearing that entry is absent and this is the plain internal energy.  Kept honest at ONE seam: the
    //! iterator's E-flat gate, the facade GetEnergy() and the display all read this (doc/GPWPlan1.md 4b).
    double GetTotalEnergy     () const;
    //! \f$\sum\f$ {Potential, Constant} -- the virial's denominator (the constants count: the theorem is for the
    //! whole system).
    double GetPotentialEnergy () const;
    //! \f$\sum\f$ {Kinetic, Potential} -- the band-structure electronic energy, EXCLUDING the structure constants
    //! (\f$E_{nn}\f$, \f$E_{\alpha Z}\f$), rest mass and entropy: for a plane-wave crystal exactly the prototype's
    //! "electronic" energy, a clean SCF stationary-point cross-check.
    double GetElectronicEnergy() const;
    //! \f$\sum\f$ {Kinetic}.
    double GetKineticEnergy   () const;
    double GetVirial          () const { return GetPotentialEnergy()/GetKineticEnergy(); }
    //! \brief THE BAND FORM: \f$E=\sum_i f_i\epsilon_i+\sum_{\rm terms}(E_{\rm term}-\mathrm{Tr}(DV_{\rm term}))\f$.
    //! The kinetic term's correction is identically 0, so this never evaluates \f$\langle T\rangle\f$ -- the
    //! point of the form.  \a sumFEps = \f$\sum_i f_i\epsilon_i\f$ from the wave function.  THROWS, naming the
    //! term, if any contribution has no \c TrDV.
    double GetBandEnergy(double sumFEps) const;

    EnergyBreakdown& operator+=(const EnergyBreakdown&);
    void Display() const;

    const std::vector<std::pair<std::string,EnergyTerm>>& Terms      () const { return itsTerms; }
    const std::vector<std::pair<std::string,double>>&     Diagnostics() const { return itsDiagnostics; }

    ChargeBreakdown charge;   //!< the charge accounting of the same density (see ChargeBreakdown)

private:
    double RoleSum(EnergyRole a) const;
    double RoleSum(EnergyRole a, EnergyRole b) const;
    std::vector<std::pair<std::string,EnergyTerm>> itsTerms;         //!< insertion-ordered; names unique
    std::vector<std::pair<std::string,double>>     itsDiagnostics;   //!< insertion-ordered; names unique
};

} //namespace
