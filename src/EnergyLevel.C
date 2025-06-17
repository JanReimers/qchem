// File: EnergyLevel.H  Energy level with degeneracy and orbital list.

#include <EnergyLevel.H>
#include <Orbitals/Orbitals.H>
#include <cassert>
#include <iostream>
#include <iomanip>

EnergyLevel::EnergyLevel(const Orbital* o)
    : e(o->GetEigenEnergy()), occ(o->GetOccupation()), degen(o->GetDegeneracy()), qns(o->GetQNs())
{};

EnergyLevel::EnergyLevel(const EnergyLevel& el) 
: e(el.e), occ(el.occ), degen(el.degen), qns(el.qns)//, orbital(el.orbital)
{}


void EnergyLevel::merge(const EnergyLevel& el)
{
    occ+=el.occ;
    // if (occ==1 && el.occ==1)
    //     std::cout << "adding first electron" << std::endl;
    // if (occ >1)
    //     std::cout << "occ=" << occ << std::endl;
    degen+=el.degen;
    assert(qns.ms==el.qns.ms);
    // Should we had lists of QNs and Orbitals?
}

void EnergyLevel::Report(std::ostream& os) const
{
    os.setf(std::ios::fixed,std::ios::floatfield);
    os << std::setw(14) << std::setprecision(8) << e 
       << " (" << std::setw(2) << std::setprecision(0) << occ 
       << "/"  << std::setw(2) << degen 
       << ") " << qns;
}

const EnergyLevel& EnergyLevels::find(const Orbital_QNs& oqns) const
{
    auto i=itsQNLevels.find(oqns);
    if (i==itsQNLevels.end())
        for (auto o:itsQNLevels)
            std::cout << o.first << std::endl;
    assert(i!=itsQNLevels.end());
    return i->second;
}

void EnergyLevels::merge(const EnergyLevels& els)
{
    for (auto& el:els) insert(el.second);
}

void EnergyLevels::merge(const EnergyLevels& els, double tol)
{
    // std::cout << "Existing levels" << std::endl;
    // Report(std::cout);
    // std::cout << "To merge levels" << std::endl;
    // els.Report(std::cout);
    
    for (auto& el:els) 
    {
        Irrep_QNs el_qns(el.second.qns);
        auto il=itsELevels.lower_bound(el.first-tol);
        auto iu=itsELevels.upper_bound(el.first+tol);
        if (il!=itsELevels.end())
        {
            Irrep_QNs il_qns(il->second.qns);
            bool symmatch =  (il_qns.SequenceIndex() == el_qns.SequenceIndex());
            // if (el.first<0.0)
            //     std::cout << std::setprecision(6) << el.first << " " << il->first << " " <<el_qns << " " << il_qns << std::endl;
            if ((!symmatch) || il==iu)
                insert(el.second);
            else
                il->second.merge(el.second);
        }
        else
            insert(el.second);
        
        
    }
    // std::cout << "Merged levels" << std::endl;
    // Report(std::cout);
}


void EnergyLevels::Report(std::ostream& os) const
{
    for (auto el:itsELevels) 
    if (el.second.occ>0)
    {
        el.second.Report(os);
        os << std::endl;
        // if (el.first>0.0) break;
    }
}

