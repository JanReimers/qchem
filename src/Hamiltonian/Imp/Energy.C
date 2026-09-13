// File: Hamiltonian/Imp/Energy.C  The energy breakdown's sums, merges and display.
module;
#include <iomanip>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>
#include <string_view>
module qchem.Energy;

namespace qchem
{

using std::cout;
using std::endl;

void EnergyBreakdown::Add(const std::string& name, double E, EnergyRole role, std::optional<double> TrDV)
{
    for (auto& [n,t] : itsTerms)
        if (n==name)
        {
            if (t.role!=role) throw std::logic_error("EnergyBreakdown::Add: '"+name+"' re-added with a different role");
            t.E+=E;
            if (t.TrDV && TrDV) *t.TrDV+=*TrDV;   // both halves know their expectation -> the merged one does
            else                t.TrDV.reset();   // either half does not -> the merged entry cannot claim to
            return;
        }
    itsTerms.emplace_back(name, EnergyTerm{E, TrDV, role});
}

void EnergyBreakdown::AddDiagnostic(const std::string& name, double v)
{
    for (auto& [n,d] : itsDiagnostics) if (n==name) { d+=v; return; }
    itsDiagnostics.emplace_back(name, v);
}

double EnergyBreakdown::operator[](std::string_view name) const
{
    for (const auto& [n,t] : itsTerms) if (n==name) return t.E;
    return 0.0;
}
double EnergyBreakdown::Diagnostic(std::string_view name) const
{
    for (const auto& [n,d] : itsDiagnostics) if (n==name) return d;
    return 0.0;
}
bool EnergyBreakdown::Has(std::string_view name) const
{
    for (const auto& [n,t] : itsTerms) if (n==name) return true;
    return false;
}

double EnergyBreakdown::RoleSum(EnergyRole a) const
{
    double s=0.0;
    for (const auto& [n,t] : itsTerms) if (t.role==a) s+=t.E;
    return s;
}
double EnergyBreakdown::RoleSum(EnergyRole a, EnergyRole b) const { return RoleSum(a)+RoleSum(b); }

double EnergyBreakdown::GetTotalEnergy() const
{
    double s=0.0;
    for (const auto& [n,t] : itsTerms) s+=t.E;
    return s;
}
double EnergyBreakdown::GetPotentialEnergy () const { return RoleSum(EnergyRole::Potential, EnergyRole::Constant); }
double EnergyBreakdown::GetElectronicEnergy() const { return RoleSum(EnergyRole::Kinetic,   EnergyRole::Potential); }
double EnergyBreakdown::GetKineticEnergy   () const { return RoleSum(EnergyRole::Kinetic); }

double EnergyBreakdown::GetBandEnergy(double sumFEps) const
{
    double corr=0.0;
    for (const auto& [n,t] : itsTerms)
    {
        if (!t.TrDV) throw std::logic_error("EnergyBreakdown::GetBandEnergy: contribution '"+n+
                                            "' has no Tr(D V) -- the band-energy form is not available for it");
        corr+=t.E-*t.TrDV;
    }
    return sumFEps+corr;
}

EnergyBreakdown& EnergyBreakdown::operator+=(const EnergyBreakdown& e1)
{
    for (const auto& [n,t] : e1.itsTerms)       Add(n, t.E, t.role, t.TrDV);
    for (const auto& [n,d] : e1.itsDiagnostics) AddDiagnostic(n, d);
    charge.lost+=e1.charge.lost;
    return *this;
}

void EnergyBreakdown::Display() const
{
    cout << endl;
    cout << "Total energy breakdown :" << endl;
    cout << "------------------------" << endl;
    for (const auto& [n,t] : itsTerms)       cout << std::left << std::setw(10) << n << ":" << t.E << endl;
    for (const auto& [n,d] : itsDiagnostics) cout << "  " << std::left << std::setw(8) << n << ":" << d << " (diagnostic)" << endl;
    cout << "------------------------" << endl << endl;
}

} //namespace
