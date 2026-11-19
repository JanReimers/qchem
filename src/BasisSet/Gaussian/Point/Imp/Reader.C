// File: BasisSet/Gaussian/Point/Imp/Reader.C  The one shell-reading loop of the Gaussian IrrepBasisSet readers, and its report rows.
module;
#include <map>
#include <memory>
#include <string>
#include <vector>
module qchem.BasisSet.Gaussian.Point.Reader;
import qchem.PeriodicTable;   // the element symbol of a row
import qchem.Reporting;       // basis.elements rows, only inside an open "basis" section

namespace qchem::BasisSet::Gaussian
{

std::vector<ReadShell> ReadShells(Reader& bsr, const Structure& cl, bool spherical)
{
    namespace rpt = qchem::report;
    const bool report=rpt::InSection("basis");
    std::map<int,size_t> atoms;                    // Z -> atoms of it (the row's count)
    if (report) for (const Site& site : cl.Sites()) ++atoms[site.Z];
    std::map<int,bool> reported;                   // Z -> its row filed (once per element, from its first site)

    std::vector<ReadShell> out;
    for (const Site& site : cl.Sites())
    {
        const size_t first=out.size();
        bsr.FindAtom(site);
        while (std::unique_ptr<GaussianRF> rf{bsr.ReadNext(site)})
        {
            const std::vector<int> newLs=bsr.GetLs();
            bool duplicate=false;
            for (auto& s : out)
                if (*s.rf==*rf && s.Ls==newLs) { duplicate=true; break; }   // the same shell twice: drop the repeat
            if (!duplicate) out.push_back({std::move(rf), newLs});          // another l on a shared radial: its own block
        }
        if (!report || reported[site.Z]) continue;
        reported[site.Z]=true;
        std::map<int,size_t> perL;                 // l -> contracted shells carrying it on this atom
        size_t nFunctions=0;
        for (size_t i=first; i<out.size(); ++i)
            for (int l : out[i].Ls)
            {
                ++perL[l];
                nFunctions+= spherical ? size_t(2*l+1) : size_t((l+1)*(l+2)/2);
            }
        std::string shells;
        for (const auto& [l,n] : perL) shells+=(shells.empty() ? "" : " ")+std::to_string(n)+(l<8 ? std::string(1,"spdfghik"[l]) : "(l="+std::to_string(l)+")");
        rpt::Row row("elements");
        rpt::Set("element",    thePeriodicTable().GetSymbol(site.Z));
        rpt::Set("Z",          (long)site.Z);
        rpt::Set("atoms",      (long)atoms[site.Z]);
        rpt::Set("shells",     shells);
        rpt::Set("nFunctions", (long)nFunctions);
    }
    return out;
}

} //namespace qchem::BasisSet::Gaussian
