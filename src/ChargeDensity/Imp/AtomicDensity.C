// File: ChargeDensity/Imp/AtomicDensity.C  SAD atomic-density database reader + radial interpolation.
module;
#include <cassert>
#include <fstream>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <vector>
#include <memory>
#include <map>
#include <nlohmann/json.hpp>

module qchem.ChargeDensity.AtomicDensity;
import qchem.Math;   // std::sqrt/log/exp/floor + Pi (qchem.CMath) + norm(Vector3D)

#ifndef CHARGEDENSITY_DATA_PATH
#error "CHARGEDENSITY_DATA_PATH must be defined by CMake"
#endif

namespace qchem::ChargeDensity
{

//! The parsed database (array of per-element-x-functional entries) for \a dbfile, loaded once on first use
//! (one cache entry per file name; construct-on-first-use, no static-init-order issue).
static const nlohmann::json& database(const std::string& dbfile)
{
    static std::map<std::string, nlohmann::json> cache;
    auto it = cache.find(dbfile);
    if (it!=cache.end()) return it->second;
    std::ifstream f(std::filesystem::path(CHARGEDENSITY_DATA_PATH) / dbfile);
    if (!f) throw std::runtime_error("AtomicDensity: cannot open " + dbfile + " under " CHARGEDENSITY_DATA_PATH);
    nlohmann::json j; f >> j;
    return cache.emplace(dbfile, std::move(j)).first->second;
}

//----------------------------------------------------------------------------------- RadialDensity

RadialDensity::RadialDensity(double rmin, double rmax, std::vector<double> rho)
    : itsMesh(qcMesh::MakeRadial({.radial=qcMesh::RadialKind::Log, .nRadial=(int)rho.size(),
                                  .logStart=rmin, .logStop=rmax}))
    , itsRho(rho.size(), rho.data())   // tabulated rho on the mesh nodes
    , itsCharge(4.0*Pi*qcMesh::Integrate(itsMesh, itsRho))   // 4*pi*int r^2 rho dr (weights fold r^2)
{
    assert(itsRho.size()==itsMesh.size());
}

double RadialDensity::operator()(double r) const
{
    const rvec_t& R=itsMesh.R();
    if (r<=R[0])            return itsRho[0];         // flat core clamp
    if (r>=R[R.size()-1])   return 0.0;              // stored tail is already ~0
    // bracket r between mesh nodes (monotonic) and linearly interpolate
    size_t lo=0, hi=R.size()-1;
    while (hi-lo>1) { size_t mid=(lo+hi)/2; if (R[mid]<=r) lo=mid; else hi=mid; }
    double f = (r-R[lo])/(R[hi]-R[lo]);
    return itsRho[lo]*(1.0-f) + itsRho[hi]*f;
}

double RadialDensity::FormFactor(double G) const
{
    // 4*pi * int rho(r) sinc(G r) r^2 dr = Sum_i 4*pi*w_i rho_i sinc(G r_i)  (w_i carries the r^2 jacobian).
    const rvec_t& R=itsMesh.R();
    const rvec_t& W=itsMesh.W();
    double sum = 0;
    for (size_t i=0;i<R.size();i++)
    {
        double x = G*R[i];
        double sinc = (x<1e-8) ? 1.0 : std::sin(x)/x;     // sinc(0)=1
        sum += W[i]*itsRho[i]*sinc;
    }
    return 4.0*Pi*sum;
}

// Match (Z, functional) and -- if Nval>=0 -- the charge state Nelec==Nval; return the entry (or end()).
static const nlohmann::json* FindAtomicEntry(int Z, const std::string& functional,
                                             const std::string& dbfile, int Nval)
{
    for (const nlohmann::json& e : database(dbfile))
        if (e.value("Z",-1)==Z && e.value("functional",std::string())==functional
            && (Nval<0 || e.value("Nelec",-1)==Nval))
            return &e;
    return nullptr;
}

RadialDensity GetAtomicDensity(int Z, const std::string& functional, const std::string& dbfile, int Nval)
{
    if (const nlohmann::json* e = FindAtomicEntry(Z, functional, dbfile, Nval))
    {
        const nlohmann::json& g = e->at("grid");
        return RadialDensity(g.at("rmin").get<double>(), g.at("rmax").get<double>(),
                             e->at("rho").get<std::vector<double>>());
    }
    throw std::runtime_error("AtomicDensity: no entry for Z=" + std::to_string(Z)
                             + " functional='" + functional + "'"
                             + (Nval>=0 ? " Nelec=" + std::to_string(Nval) : std::string())
                             + " in " + dbfile);
}

bool HasAtomicDensity(int Z, const std::string& functional, const std::string& dbfile, int Nval)
{
    return FindAtomicEntry(Z, functional, dbfile, Nval) != nullptr;
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
