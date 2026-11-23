module;
#include <ostream>   // R1.9: free operator<<(ostream&,const char*) is NOT visible via a module import
#include <iomanip>
#include <cassert>
#include <memory>
#include <vector>
#include <stdexcept>
#include <string>

module qchem.Structure;
import qchem.Streamable;
import qchem.Symmetry.Molecule.AbelianGroup;   // BuildAbelianGroup, Centroid (MakeSymmetry)

namespace qchem {

Molecule::Molecule(const Structure& atoms)
{
    for (const Site& a:atoms.Sites()) itsSites.push_back(a);
}

void Molecule::Insert(const Site& a)
{
    assert(a.Z>0 && a.Z<150);   // the Atom ctor's range check, now that a Site needs no Atom to be stored
    itsSites.push_back(a);
}

std::ostream& Molecule::Write(std::ostream& os) const
{
    os << "Molecule with " << GetNumAtoms() << " atoms"
    << ", nuclear charge " << GetNuclearCharge() << "(e)"
    << ", net charge "<< GetNetCharge() << "(e)" << std::endl;
    os << "Atom #  Element  Position vector" << std::endl;
    int i=1;
    for (const Site& b:Sites()) {os << std::setw(5) << i++ << "   "; WriteSite(os,b) << std::endl;}
    // os << std::endl;

    return os;
}


// Moved from the interface (D-BUILD-BLAST, lane B 2026-11-18): an edit here recompiles this TU, not the ~450 importers of qchem.Structure.
Molecule::Molecule() {}
std::shared_ptr<Structure> Molecule::Clone() const {return std::make_shared<Molecule>(*this);}
size_t Molecule::GetNumAtoms() const {return itsSites.size();}
Symmetry::GroupPtr Molecule::MakeSymmetry(const Symmetry::GroupTolerance& tol) const
{
    namespace SM=Symmetry::Molecule;
    if (itsSites.empty()) throw std::invalid_argument("Molecule::MakeSymmetry: an empty site set -- its group would be every isometry, "
                                                      "and there is no finite group to issue");
    std::vector<SM::SymPoint> pts;
    for (const Site& s:itsSites) pts.push_back({s.Z, s.R});
    // The ABELIAN subgroup (see the declaration): BuildAbelianGroup's ops act about the centroid of the points it was given.  The group
    // NAMES itself as that subgroup of the detected full group (user, 2026-11-20), so a report never passes it off as the point group.
    const SM::PointGroup pg=SM::DetectPointGroup(pts, tol.length);
    const std::string name=pg.abelian==pg.symbol ? pg.symbol+" (the detected point group; abelian)"
                                                 : pg.abelian+" (abelian subgroup of "+pg.symbol+")";
    return Symmetry::FiniteGroup(SM::BuildAbelianGroup(pts, tol.length).ops, SM::Centroid(pts), tol, name);
}
const Site& Molecule::GetSite(size_t i) const {return itsSites[i];}
      Site& Molecule::GetSite(size_t i)       {return itsSites[i];}

} // namespace qchem