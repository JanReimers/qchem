// File: EnergyLevel.C  Energy level with degeneracy and orbital list.
module;
#include <cassert>
#include <iostream>
#include <iomanip>
module qchem.EnergyLevel;
import qchem.Streamable;

namespace qchem::Orbitals
{

EnergyLevel::EnergyLevel(const Orbital* o)
    : e(o->GetEigenEnergy()), occ(o->GetOccupation()), degen(o->GetDegeneracy()), qns(o->GetQNs())
{
    // REPORTING scale (Symmetry::StarSize doc): an IBZ wedge block stands for its whole k-star, so its
    // levels report full-mesh occ/degen (e.g. the X-star shell as 8/8, matching the unfolded run's merged
    // table) -- 1 for everything but a folded Bloch block.  Filling/weights upstream are untouched.
    const size_t star=qns.sym->StarSize();
    if (star>1) { occ*=double(star); degen*=int(star); }
};

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
    {
        std::cout << "EnergyLevels cannot find:" << oqns << " int this list:" << std::endl;
        for (auto o:itsQNLevels)
            std::cout << "   " << o.first << std::endl;
    }
    assert(i!=itsQNLevels.end());
    return i->second;
}

void EnergyLevels::merge(const EnergyLevels& els)
{
    for (auto& el:els) insert(el.second);
}

void EnergyLevels::merge(const EnergyLevels& els, double tol)
{
    for (auto& el:els)
    {
        const Irrep el_qns(el.second.qns);
        // Scan the WHOLE ±tol window for a mergeable host (the old code tested only the FIRST level in
        // the window -- an order-dependent miss).  A host is: the SAME irrep (the historical rule), or --
        // when BOTH symmetries opt in (Symmetry::MergeAcrossIrreps, i.e. crystal k-points) -- any level of
        // the same spin channel: a k-STAR of band levels is one degenerate shell of the crystal group, and
        // MergeTol's contract is to fuse it (display/cfg layer; the occupation decisions are per-block,
        // upstream in tIrrepWF::FillOrbitals).
        auto il=itsELevels.lower_bound(el.first-tol);
        auto iu=itsELevels.upper_bound(el.first+tol);
        auto host=itsELevels.end();
        for (auto i=il; i!=iu; ++i)
        {
            const Irrep i_qns(i->second.qns);
            const bool symmatch  = (i_qns.SequenceIndex() == el_qns.SequenceIndex());
            const bool starmatch = el_qns.sym->MergeAcrossIrreps() && i_qns.sym->MergeAcrossIrreps()
                                && i->second.qns.ms == el.second.qns.ms;
            if (symmatch || starmatch) { host=i; break; }
        }
        if (host!=itsELevels.end()) host->second.merge(el.second);
        else                        insert(el.second);
    }
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

} //namespace