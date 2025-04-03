// File: EnergyLevel.H  Energy level with degeneracy and orbital list.

#include <EnergyLevel.H>
#include <Symmetry.H>
#include <cassert>
#include <iostream>
#include <iomanip>

void EnergyLevel::merge(const EnergyLevel& el)
{
    occ+=el.occ;
    // if (occ==1 && el.occ==1)
    //     std::cout << "adding first electron" << std::endl;
    // if (occ >1)
    //     std::cout << "occ=" << occ << std::endl;
    degen+=el.degen;
    assert(s==el.s);
    // Should we had lists of QNs and Orbitals?
}

void EnergyLevel::Report(std::ostream& os) const
{
    os.setf(std::ios::fixed,std::ios::floatfield);
    os << std::setw(14) << std::setprecision(8) << e 
       << " (" << std::setw(2) << std::setprecision(0) << occ 
       << "/"  << std::setw(2) << degen 
       << ") " << *qn;
}


void EnergyLevels::merge(const EnergyLevels& els)
{
    for (auto& el:els) itsELevels.insert(el);
}

void EnergyLevels::merge(const EnergyLevels& els, double tol)
{
    for (auto& el:els) 
    {
        auto il=itsELevels.lower_bound(el.first-tol);
        auto iu=itsELevels.upper_bound(el.first+tol);
        if (il==itsELevels.end() || il==iu)
            itsELevels.insert(el);
        else
            il->second.merge(el.second);
        
    }
}


void EnergyLevels::Report(std::ostream& os) const
{
    for (auto el:itsELevels) 
    {
        el.second.Report(os);
        os << std::endl;
        if (el.first>0.0) break;
    }
}

