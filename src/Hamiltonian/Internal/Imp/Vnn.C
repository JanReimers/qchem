// File: Vnn.C  Nuclear-Nuclear potential.
module;
#include <iostream>
#include <cassert>
#include <vector>

module qchem.Hamiltonian.Internal.Terms;
import qchem.Energy;
import qchem.Structure;
import qchem.UnitCell;   // periodic ion-ion -> Ewald lattice sum
import qchem.Ewald;      // EwaldEnergy
import qchem.Types;      // rvec_t
import qchem.Blaze;

namespace qchem::Hamiltonian
{

Vnn::Vnn(const st_t& st)
    : Static_HT_Imp()
    , theStructure(st)
{};

rsmat_t Vnn::CalculateMatrix(const obs_t* bs,const Spin&) const
{
    int n=bs->GetNumFunctions();
    rsmat_t ret=blazem::zero<double>(n);
    return ret;
}

void Vnn::GetEnergy(EnergyBreakdown& te,const DM_CD* cd) const
{
    if (!theStructure->isFinite())
    {
        // Periodic cell: the bare pair sum is only conditionally convergent, so use the Ewald lattice
        // sum.  The ion charges are the atoms' itsZ (for a pseudopotential crystal the cell carries the
        // VALENCE/ion charge Zion there, the natural charge of the ion-ion term).
        const UnitCell* cell=dynamic_cast<const UnitCell*>(theStructure.get());
        assert(cell && "Vnn: a non-finite Structure must be a UnitCell for the Ewald lattice sum");
        rvec_t q(cell->GetNumAtoms());
        for (size_t a=0; a<cell->GetNumAtoms(); a++) q[a]=(*cell)[a]->itsZ;
        te.Enn=EwaldEnergy(*cell,q);
        return;
    }

    double vnn=0.0;
    for(const auto& atom1:*theStructure)
        for(const auto& atom2:*theStructure)
        {
            Vector3D<double> r1=atom1->itsR, r2=atom2->itsR;
            if (r1!=r2) vnn += 0.5 * atom1->itsZ * atom2->itsZ / norm(r1-r2);
        }

    te.Enn=vnn;
}

std::ostream& Vnn::Write(std::ostream& os) const
{
    size_t Na=theStructure->GetNumAtoms();
    if (theStructure->isFinite())
        os << "    Nuclear-Nuclear potential ZiZj/|Ri-Rj| with " << Na*(Na-1) << " nucleus pairs." << std::endl;
    else
        os << "    Ion-Ion (Ewald lattice sum) over " << Na << " ions per cell." << std::endl;
    return os;
}

} //namespace
