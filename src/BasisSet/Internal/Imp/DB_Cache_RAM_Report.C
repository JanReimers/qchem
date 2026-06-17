// File: BasisSet/Internal/Imp/DB_Cache_RAM.C
module;
#include <iostream>
#include <iomanip>
#include <cassert>
#include <map>
#include <chrono>
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

template <class T>  void IntegralsCache_RAM<T>::ReportRAMUsage(std::ostream& os) const
{
    bool verbose=false;
    os << "------------------- IntegralsCache RAM usage report -------------------------" << endl;
    size_t J_ram=Report(Jac,"Direct   Jac",verbose); //# of Ts
    size_t K_ram=Report(Kab,"Exchange Kab",verbose);

    size_t cach4_ram=0;
    for (auto& i:itsCache4s)
       cach4_ram=i.second->RAMsize();
    size_t total=J_ram+K_ram+cach4_ram;

    os << "  Jac   data: "  << ram(J_ram,sizeof(T)) << percent(J_ram,total) << endl;
    os << "  Kab   data: "  << ram(K_ram,sizeof(T)) << percent(K_ram,total) << endl;
    os << "  Cach4 data: "  << ram(cach4_ram,sizeof(T)) << percent(cach4_ram,total) << endl;
}

std::ostream& operator << (std::ostream& os, const IntegralsCache_Base::IBS_ID_t& ids)
{
    return os << "radial=" << std::get<0>(ids) << " angular=" << std::get<1>(ids);
}

std::ostream& operator << (std::ostream& os, const std::pair<IntegralsCache_Base::IBS_ID_t,IntegralsCache_Base::IBS_ID_t>& ids)
{
    return os << " a=" << ids.first << "   b=" << ids.second;
}
template <class T>  size_t IntegralsCache_RAM<T>::Purge(map4_t& eri4s,const id_pair_t& old)
{
    size_t ret=0;
    auto i1=eri4s.find(old.first);
    if (i1!=eri4s.end())
    {
        auto i2=i1->second.find(old.second);
        if (i2!=i1->second.end())
        {
            bool ab1=i1->first ==itsLastKey4a;
            bool ab2=i2->first==itsLastKey4b;
            bool ba1=i1->first ==itsLastKey4b;
            bool ba2=i2->first==itsLastKey4a;
            // itsLogger << "bools=" << ab1 << " " << ab2 << " " << ba1 << " " << ba2 << endl;
            assert(!(ab1&&ab2));
            assert(!(ba1&&ba2));
            ret+=i2->second.size()*sizeof(T)/1024/1024;
            i1->second.erase(i2);
            if (i1->second.empty()) eri4s.erase(i1);
            itsLogger << "GC purging " << ret << "(MB) for "<< old << endl;
        }
        
    }
    return ret;
}
template <class T>  void IntegralsCache_RAM<T>::RunGarbageCollector()
{
    // itsLogger << "Running IntegralsCache_RAM garbage collector:" << endl;
    // ReportRAMUsage(itsLogger);
    // itsLogger << "  Total ERI4 RAM=" << itsTotalRAM << "(MB)" << endl;

    size_t N=0;
    while (itsTotalRAM>itsMaxRAM && N<100)
    {
        // itsLogger << "  ERI4_timestamps list:" << endl;
        // for (auto& i:ERI4_timestamps)
        // {
        //     itsLogger << i.first << " " << i.second << endl;
        // }
        std::map<time_t,id_pair_t> oldest;
        for (auto& i:ERI4_timestamps) oldest[i.second]=i.first;
        // itsLogger << "  Oldest first list:" << endl;
        // for (auto& i:oldest)
        // {
        //     itsLogger << i.first << " " << i.second << endl;
        // }


        id_pair_t oldab = oldest.begin()->second;
        // itsLogger << "old=" << oldab;
        id_pair_t oldba=std::make_pair(oldab.second,oldab.first);
        size_t ram=0;
        bool ab1=oldab.first==itsLastKey4a;
        bool ab2=oldab.second==itsLastKey4b;
        bool ba1=oldab.first==itsLastKey4b;
        bool ba2=oldab.second==itsLastKey4a;
        // itsLogger << "bools=" << ab1 << " " << ab2 << " " << ba1 << " " << ba2 << endl;
        if ((ab1&&ab2) || (ba1&&ba2))
        {
            //  itsLogger << "   Purge itsLastKey4 a&&b clash." << endl;
        }
        else
        {
            ram+=Purge(Jac,oldab);
            ram+=Purge(Kab,oldab);
            ram+=Purge(Jac,oldba);
            ram+=Purge(Kab,oldba);
            assert(ram>0);
            {   
                ERI4_timestamps.erase(oldab);
                ERI4_timestamps.erase(oldba);
                assert(ram<itsTotalRAM);
                itsTotalRAM-=ram;
                // itsLogger << "  Total ERI4 RAM=" << itsTotalRAM << "(MB) " << oldab << endl;
            }
        }
        N++;
        
    }
    // ReportRAMUsage(itsLogger);
}



template struct IntegralsCache_RAM<double>;
} //namespace
