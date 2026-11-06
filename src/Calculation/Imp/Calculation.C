// File: Calculation/Imp/Calculation.C  What every facade shares.
module;
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
module qchem.Calculation;
import qchem.Reporting;

namespace qchem
{

void Calculation::ReportResult(const std::string& label, bool converged, size_t iterations, double charge, const EnergyBreakdown& E,
                               double energyTarget)
{
    namespace rpt = qchem::report;
    const double dq=E.Diagnostic("HartreeChargeError"), dE=E.Diagnostic("HartreeChargeEnergy");
    const bool gridHartree=dq!=0.0;   // only a grid route writes it, and a mesh never integrates a density exactly
    {
        rpt::Section result("result");
        rpt::Set("converged",  converged);
        rpt::Set("iterations", (long)iterations);
        rpt::Set("charge",     charge);
        rpt::Set("Eelec",      E.GetElectronicEnergy());
        rpt::Set("Etot",       E.GetTotalEnergy());
        if (gridHartree)
        {
            rpt::Set("hartreeChargeError",  dq);
            rpt::Set("hartreeChargeEnergy", dE);
            if (energyTarget>0.0) rpt::Set("hartreeChargeEnergy/target", std::fabs(dE)/energyTarget);
        }
    }
    // Precision STATED, not inherited from whatever the trace left on cout: doc/Benchmark.md compares codes at the 1e-5 Ha level,
    // and anchordiff reads Etot at 1e-10.
    const std::ios_base::fmtflags flags=std::cout.flags();
    const std::streamsize         prec =std::cout.precision();
    std::cout<<std::fixed<<std::setprecision(11)
             <<"["<<label<<"] "<<(converged ? "CONVERGED" : "NOT converged")<<" iters="<<iterations<<" charge="<<charge
             <<" Eelec="<<E.GetElectronicEnergy()<<" Etot="<<std::setprecision(10)<<E.GetTotalEnergy()<<std::setprecision(11)<<"  (";
    const char* sep="";
    for (const auto& [name,shown] : {std::pair{"Kinetic","Ekin"}, {"Een","Een"}, {"Eee","Eee"}, {"Exc","Exc"}, {"Enn","Enn"}, {"E_alphaZ","E_alphaZ"}})
        if (E.Has(name)) { std::cout<<sep<<shown<<"="<<E[name]; sep=" "; }
    std::cout<<")"<<std::endl;
    if (gridHartree)
    {
        std::cout<<std::scientific<<std::setprecision(3)<<"["<<label<<"] Hartree mesh charge error="<<dq<<" energy="<<dE;
        if (energyTarget>0.0) std::cout<<" target="<<energyTarget<<" ratio="<<std::fabs(dE)/energyTarget;
        std::cout<<std::endl;
    }
    std::cout.flags(flags);
    std::cout.precision(prec);
}

} // namespace qchem
