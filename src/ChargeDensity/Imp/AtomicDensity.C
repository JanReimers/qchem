// File: ChargeDensity/Imp/AtomicDensity.C  SAD atomic-density database reader + radial interpolation.
module;
#include <fstream>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <vector>
#include <memory>
#include <nlohmann/json.hpp>

module qchem.ChargeDensity.AtomicDensity;
import qchem.Math;   // std::sqrt/log/exp/floor (qchem.CMath) + norm(Vector3D)

#ifndef CHARGEDENSITY_DATA_PATH
#error "CHARGEDENSITY_DATA_PATH must be defined by CMake"
#endif

namespace qchem::ChargeDensity
{

//! The parsed database (array of per-element-x-functional entries), loaded once on first use.
static const nlohmann::json& database()
{
    static const nlohmann::json db = []
    {
        std::ifstream f(std::filesystem::path(CHARGEDENSITY_DATA_PATH) / "atomic_densities.json");
        if (!f) throw std::runtime_error("AtomicDensity: cannot open atomic_densities.json under "
                                         CHARGEDENSITY_DATA_PATH);
        nlohmann::json j; f >> j; return j;
    }();
    return db;
}

//----------------------------------------------------------------------------------- RadialDensity

RadialDensity::RadialDensity(double rmin, double rmax, std::vector<double> rho)
    : itsRmin(rmin), itsRmax(rmax)
    , itsLogStep(std::log(rmax/rmin)/(rho.size()-1))
    , itsRho(std::move(rho))
    , itsCharge(0.0)
{
    // 4*pi*int r^2 rho dr on the log grid: dr = r du, du = itsLogStep, trapezoid (half weight at the ends).
    const int N = (int)itsRho.size();
    for (int i=0;i<N;i++)
    {
        double r = itsRmin*std::exp(i*itsLogStep);
        double w = (i==0||i==N-1) ? 0.5 : 1.0;
        itsCharge += w * 4.0*Pi*r*r*itsRho[i] * r*itsLogStep;
    }
}

double RadialDensity::operator()(double r) const
{
    if (r<=itsRmin) return itsRho.front();   // flat core clamp
    if (r>=itsRmax) return 0.0;              // stored tail is already ~0
    double u = std::log(r/itsRmin)/itsLogStep;
    int    i = (int)std::floor(u);
    if (i>=(int)itsRho.size()-1) return itsRho.back();
    double f = u-i;
    return itsRho[i]*(1.0-f) + itsRho[i+1]*f;
}

RadialDensity GetAtomicDensity(int Z, const std::string& functional)
{
    for (const nlohmann::json& e : database())
        if (e.value("Z",-1)==Z && e.value("functional",std::string())==functional)
        {
            const nlohmann::json& g = e.at("grid");
            return RadialDensity(g.at("rmin").get<double>(), g.at("rmax").get<double>(),
                                 e.at("rho").get<std::vector<double>>());
        }
    throw std::runtime_error("AtomicDensity: no entry for Z=" + std::to_string(Z)
                             + " functional='" + functional + "' in atomic_densities.json");
}

//-------------------------------------------------------------------------------- RecentredAtomicDensity

RecentredAtomicDensity::RecentredAtomicDensity(std::shared_ptr<const RadialDensity> rad, const rvec3_t& R)
    : itsRadial(std::move(rad)), itsR(R)
{
}

double RecentredAtomicDensity::operator()(const rvec3_t& r) const
{
    return (*itsRadial)(norm(r-itsR));
}

rvec3_t RecentredAtomicDensity::Gradient(const rvec3_t& r) const
{
    const double h=1e-4;
    auto d=[&](const rvec3_t& a, const rvec3_t& b){ return ((*this)(a)-(*this)(b))/(2.0*h); };
    return rvec3_t(d(rvec3_t(r.x+h,r.y,r.z), rvec3_t(r.x-h,r.y,r.z)),
                   d(rvec3_t(r.x,r.y+h,r.z), rvec3_t(r.x,r.y-h,r.z)),
                   d(rvec3_t(r.x,r.y,r.z+h), rvec3_t(r.x,r.y,r.z-h)));
}

} //namespace
