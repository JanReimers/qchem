// File: BasisSet1/Internal/Imp/DB_Cache_RAM.C
module;
#include <iostream>
#include <iomanip>
module qchem.BasisSet.Internal.DB_Cache_RAM;

namespace BasisSet
{

using std::cout;
using std::endl;

template <class T> size_t IntegralsCache_RAM<T>::Report(const map4_t& m4, const std::string& name, bool verbose)
{
    size_t nt=0; //# of Ts stored.  *szie_of(T) latter.
    if (verbose) cout << "Report for " << name << ":" << endl;
    for (auto i:m4)
    for (auto j:i.second)
    {
        nt+=j.second.size();
        if (verbose)
            cout << "        "  << j.second.size() << "  i=" << std::get<1>(i.first) << "  j=" << std::get<1>(j.first) << endl;
    }
    if (verbose) cout << "  Total:" << nt << endl;
    return nt;
}

std::string percent(size_t n, size_t total)
{
    double p= (100.0*n)/total;
    std::ostringstream os;
    os << " " << std::setprecision(0) << std::fixed << std::setw(3) << p << "%";
    return os.str();
}

std::string ram(size_t n, size_t tsize)
{
    size_t megas= n*tsize/1024/1024;
    std::ostringstream os;
    os << " " << std::setprecision(0) << std::fixed << std::setw(12) << megas << " (MB)";
    return os.str();
}

template <class T>  void IntegralsCache_RAM<T>::ReportRAMUsage() const
{
    bool verbose=false;
    cout << "------------------- IntegralsCache RAM usage report -------------------------" << endl;
    size_t J_ram=Report(Jac,"Direct   Jac",verbose);
    size_t K_ram=Report(Kab,"Exchange Kab",verbose);

    size_t cach4_ram=0;
    for (auto& i:itsCache4s)
       cach4_ram=i.second->RAMsize();
    size_t total=J_ram+K_ram+cach4_ram;

    cout << "  Jac   data: "  << ram(J_ram,sizeof(T)) << percent(J_ram,total) << endl;
    cout << "  Kab   data: "  << ram(K_ram,sizeof(T)) << percent(K_ram,total) << endl;
    cout << "  Cach4 data: "  << ram(cach4_ram,sizeof(T)) << percent(cach4_ram,total) << endl;
}

template struct IntegralsCache_RAM<double>;
} //namespace
