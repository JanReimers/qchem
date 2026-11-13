// File: ChargeDensity/AtomicDensity.C  Reader for the SAD atomic-density database (radial rho(r) per element).
//
// The database (Data/atomic_densities.json) is generated offline by `scfrun --model LDA --out ...`: one
// neutral-atom radial density per element x functional, sampled on a log-radial grid (see
// doc/SCFSeedingPlan.md section 9).  This module reads it and exposes the radial density two ways: as a 1-D
// RadialDensity (interpolated rho(r)), and -- recentred at a nucleus -- as a 3-D ScalarFunction<double> the
// scalar function-fitter can sample to build the molecular SAD seed.
module;
#include <string>
#include <vector>
#include <memory>
#include <utility>
#include <map>
export module qchem.ChargeDensity.AtomicDensity;
import qchem.ScalarFunction;   // ScalarFunction<double>, rvec3_t
import qchem.Mesh.Quadrature;  // qcMesh::RadialMesh, MakeRadial, Integrate (radial quadrature)

export namespace qchem::ChargeDensity
{

//! WHICH seed densities a run asks the library for: the XC \a functional and, per atomic number, the
//! pseudopotential variant \c q (valence electron count of the GTH potential the run's Hamiltonian uses).
//! The run setup builds it from its \c PseudopotentialSet, so the seed, the terms and the basis cannot
//! disagree about which pseudo-atom they are (C4 step 3; Li q1 and q3 are different seed densities).  A Z with no
//! \c q entry is unkeyed (today's neutral-discovery rule, which THROWS when two variants are neutral).
//! Implicitly built from a bare functional name, which keeps unkeyed callers (tests, molecular) unchanged.
struct SeedKey
{
    std::string      functional = "LDA";
    std::map<int,int> q;                                   //!< Z -> pseudopotential valence count
    //! The IonicSAD seed's per-SITE formal charges, in atom order, when the run STATES them (\c SolidCalcOptions::siteCharges); EMPTY =
    //! derived by \c IonicFormalCharges, whose tie-break puts the charge states of symmetry-equivalent sites by ATOM INDEX (a mixed-
    //! valence spinel's Mn3+/Mn4+ order is then a choice nobody made).
    std::vector<int>  siteCharges;
    //! The polarized seed's per-SITE spins (+1, -1, 0), in atom order, when the run STATES its magnetic ordering
    //! (\c SolidCalcOptions::siteSpins); EMPTY = the structure's per-atom flip bits.  ONE source with the imposition
    //! (OpenWork §3 "One magnetic ordering has TWO sources of truth", 2026-11-12): an inline cell without atom spins
    //! and a stated AFM ordering once seeded ferromagnetic, which the imposed Shubnikov group then projected to m = 0.
    std::vector<int>  siteSpins;
    SeedKey() = default;
    SeedKey(const char* f) : functional(f) {}
    SeedKey(std::string f, std::map<int,int> qByZ={}) : functional(std::move(f)), q(std::move(qByZ)) {}
    int Q(int Z) const {auto it=q.find(Z); return it==q.end() ? -1 : it->second;}   //!< -1 = unkeyed
};

//! THE SEED LIBRARY a run reads for \a Z, one rule for every structure (OpenWork §3 "A pseudopotential run seeds from the VALENCE
//! library for its q, on every structure", lane B 2026-11-12): a Z the key gives a pseudopotential \c q is the PSEUDO-atom the
//! Hamiltonian's terms use, so it reads the VALENCE library at that q (the core the pseudopotential removed is not in the run); an
//! unkeyed Z (an all-electron run) reads the ALL-ELECTRON library.  Before this a PP molecule seeded from the all-electron densities,
//! core charge included, which the fit's charge constraint then squeezed to the valence count.
inline std::string SeedLibrary(const SeedKey& key, int Z)
{
    return key.Q(Z)>=0 ? "atomic_valence_densities.json" : "atomic_densities.json";
}

//! A neutral atom's spherically-averaged radial density rho(r), read from the SAD database.  It lives on a
//! qcMesh::RadialMesh (a Log mesh: nodes r_i = rmin*q^i, with the r^2 jacobian folded into the weights), so
//! the integrals below are plain weighted sums.  operator()(r) linearly interpolates between nodes; outside
//! [rmin,rmax] it clamps (rho(rmin) for the core, 0 beyond rmax -- the stored tail is already ~0).
class RadialDensity
{
public:
    RadialDensity(double rmin, double rmax, std::vector<double> rho);   //!< builds a Log RadialMesh internally
    double operator()(double r) const;   //!< interpolated rho(r)
    double Charge() const {return itsCharge;}   //!< 4*pi*int r^2 rho dr (= Sum_i 4*pi*w_i rho_i; ~ Nelec)
    int    GetN() const {return (int)itsMesh.size();}
    //! The table's SUPPORT radius: rho(r)==0 exactly for r>=Rmax() (the clamp in operator()).  A periodic
    //! consumer's lattice-image sum is therefore FINITE and exact -- it needs only the images whose support
    //! ball reaches the point (SeedCD's periodic rho(r)); no distance cut is involved.
    double Rmax() const {return itsMesh.R()[itsMesh.size()-1];}
    //! Reciprocal-space form factor \f$\tilde\rho(G)=4\pi\int_0^\infty \rho(r)\,\mathrm{sinc}(Gr)\,r^2\,dr\f$
    //! (the radial Fourier transform of a spherical density).  \f$\tilde\rho(0)=\f$ Charge().  Used by the
    //! plane-wave SAD seed's structure-factor sum.
    double FormFactor(double G) const;
private:
    qcMesh::RadialMesh  itsMesh;     //!< Log radial mesh: nodes R()[i] + r^2-folded weights W()[i]
    rvec_t              itsRho;      //!< rho at each mesh node
    double              itsCharge;   //!< cached 4*pi*int r^2 rho dr
};

//! Read element \a Z's radial density for \a functional from database \a dbfile (throws if absent).  The
//! default holds all-electron densities (the molecular SAD source); pass "atomic_valence_densities.json"
//! for the pseudo-valence densities (the plane-wave SAD source).  \a Nval selects a CHARGE STATE by valence
//! electron count (the entry's \c Nelec): \f$<0\f$ (default) = the first \a Z match (the neutral atom);
//! \f$\ge0\f$ = the entry with \c Nelec==Nval (e.g. \c Nval=8 for F\f$^-\f$).  The seed-density library holds
//! \a q >= 0 keeps only the entries generated under PP variant q (the entry's \c q), which resolves the Li q1/q3 ambiguity.
//! The library holds neutral + chemically-plausible ion entries (generated offline by \c qchem::ValenceBasisGen).
RadialDensity GetAtomicDensity(int Z, const std::string& functional="LDA",
                               const std::string& dbfile="atomic_densities.json", int Nval=-1, int q=-1);
//! Is there an entry for \a Z / \a functional (and \c Nelec==Nval if \a Nval>=0) in \a dbfile?  Lets a caller
//! (IonicSAD) prefer a charge-state density when the library has it and fall back cleanly when it does not.
bool          HasAtomicDensity(int Z, const std::string& functional="LDA",
                               const std::string& dbfile="atomic_densities.json", int Nval=-1, int q=-1);

//! SPIN-RESOLVED pair {rho_up, rho_dn} for an open-shell (Hund ground state) entry -- the spin-polarized
//! SAD source (doc/SCFSeedingPlan.md section 10).  Stored in the UP-MAJORITY convention: \c first is the
//! MAJORITY channel (\f$\int\rho_\uparrow \ge \int\rho_\downarrow\f$); the magnetic CONFIGURATION (which
//! sites flip) is assembly-time data, never in the tables.  The entry's \c rho remains the spin SUM, so
//! spin-agnostic readers (GetAtomicDensity) are untouched.  Entry selection is identical to
//! GetAtomicDensity; throws if the entry is absent or carries no pair (closed shells store none) --
//! callers gate on HasAtomicSpinPair and fall back to the rho/2 split.
std::pair<RadialDensity,RadialDensity>
              GetAtomicSpinPair(int Z, const std::string& functional="LDA",
                                const std::string& dbfile="atomic_valence_densities.json", int Nval=-1, int q=-1);
//! Does the (Z, functional[, Nelec==Nval]) entry exist AND carry a spin-resolved pair?
bool          HasAtomicSpinPair(int Z, const std::string& functional="LDA",
                                const std::string& dbfile="atomic_valence_densities.json", int Nval=-1, int q=-1);

//! A RadialDensity recentred at a nucleus \a R: a 3-D ScalarFunction rho(|r-R|) for the function-fitter.
class RecentredAtomicDensity : public virtual ScalarFunction<double>
{
public:
    RecentredAtomicDensity(std::shared_ptr<const RadialDensity>, const rvec3_t& R);
    virtual double  operator()(const rvec3_t&) const;
    virtual rvec3_t Gradient  (const rvec3_t&) const;   //!< numerical (central difference)
    const rvec3_t& Center() const {return itsR;}              //!< the nucleus this density is centred on
    double         Rmax  () const {return itsRadial->Rmax();} //!< support radius (see RadialDensity::Rmax)
private:
    std::shared_ptr<const RadialDensity> itsRadial;
    rvec3_t                              itsR;
};

} //namespace
