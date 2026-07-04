// File: Vnn.C  Nuclear-Nuclear potential.
module;
#include <iostream>
#include <cassert>
#include <vector>
#include <functional>
#include <utility>

module qchem.Hamiltonian.Internal.Terms;
import qchem.Energy;
import qchem.Structure;
import qchem.Ewald;      // NuclearRepulsion (pair sum for finite, Ewald for periodic)
import qchem.Blaze;

namespace qchem::Hamiltonian
{

Vnn::Vnn(const st_t& st)
    : Vnn(st, [](int Z){return double(Z);})   // all-electron default: ion charge IS the true nuclear Z
{};

Vnn::Vnn(const st_t& st, std::function<double(int)> zionOf)
    : rStatic_HT_Imp()
    , theStructure(st)
    , itsZionOf(std::move(zionOf))
{
    assert(itsZionOf && "Vnn: a Z->ion-charge map is required");
};

rsmat_t Vnn::CalculateMatrix(const robs_t* bs,const Spin&) const
{
    int n=bs->GetNumFunctions();
    rsmat_t ret=blazem::zero<double>(n);
    return ret;
}

void Vnn::GetEnergy(EnergyBreakdown& te,const rDM_CD* cd) const
{
    // Pair sum for a finite molecule, Ewald lattice sum for a periodic cell (chosen by isFinite()).
    // Charges come from the Z->ion map: itsZ for all-electron, Zion for a pseudopotential.
    te.Enn=NuclearRepulsion(*theStructure, itsZionOf);
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
