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
export module qchem.ChargeDensity.AtomicDensity;
import qchem.ScalarFunction;   // ScalarFunction<double>, rvec3_t
import qchem.Mesh.Quadrature;  // qcMesh::RadialMesh, MakeRadial, Integrate (radial quadrature)

export namespace qchem::ChargeDensity
{

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
//! for the pseudo-valence densities (the plane-wave SAD source).
RadialDensity GetAtomicDensity(int Z, const std::string& functional="LDA",
                               const std::string& dbfile="atomic_densities.json");

//! A RadialDensity recentred at a nucleus \a R: a 3-D ScalarFunction rho(|r-R|) for the function-fitter.
class RecentredAtomicDensity : public virtual ScalarFunction<double>
{
public:
    RecentredAtomicDensity(std::shared_ptr<const RadialDensity>, const rvec3_t& R);
    virtual double  operator()(const rvec3_t&) const;
    virtual rvec3_t Gradient  (const rvec3_t&) const;   //!< numerical (central difference)
private:
    std::shared_ptr<const RadialDensity> itsRadial;
    rvec3_t                              itsR;
};

} //namespace
