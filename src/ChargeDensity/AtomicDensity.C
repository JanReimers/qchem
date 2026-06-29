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

export namespace qchem::ChargeDensity
{

//! A neutral atom's spherically-averaged radial density rho(r), read from the SAD database and stored on a
//! log-radial grid.  operator()(r) linearly interpolates; outside [rmin,rmax] it clamps (rho(rmin) for the
//! core, 0 beyond rmax -- the stored tail is already ~0).
class RadialDensity
{
public:
    RadialDensity(double rmin, double rmax, std::vector<double> rho);
    double operator()(double r) const;   //!< interpolated rho(r)
    double Charge() const {return itsCharge;}   //!< 4*pi*int r^2 rho dr from the stored grid (~ Nelec)
    int    GetN() const {return (int)itsRho.size();}
private:
    double              itsRmin, itsRmax, itsLogStep;   //!< log grid: r_i = rmin*exp(i*itsLogStep)
    std::vector<double> itsRho;
    double              itsCharge;
};

//! Read element \a Z's radial density for \a functional from the database (throws if absent).
RadialDensity GetAtomicDensity(int Z, const std::string& functional="LDA");

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
