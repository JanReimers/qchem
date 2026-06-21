// File: BasisSet/Molecule/PolarizedGaussian/Imp/Orbital_IBS.C  Polarized Gaussian 2-centre fit integrals.
module;
#include <cassert>

module qchem.BasisSet.Molecule.PolarizedGaussian;
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.GaussianRF;
import qchem.BasisSet.Internal.IntegralEnums;
import qchem.Cluster;
import qchem.Blaze;

// The orbital 1E / 3C (DFT) / 4C (HF) integral building now lives entirely in the Molecule-generic,
// evaluator-templated mixins (qchem.BasisSet.Molecule.IBS), instantiated with PG_Evaluator.  Nothing
// PG-specific remains here for the orbital path.  Only the EFit_IBS 2-centre fit integrals stay below:
// they use the raw radials directly (Overlap2C / Repulsion2C), not the 1E/3C/4C evaluator concepts.
namespace BasisSet::Molecule::PolarizedGaussian
{

rsmat_t MakeIntegrals(IType t2C,const PGData* ab, const Cluster* cl)
{
    assert(ab);
    int N=ab->size();
    rsmat_t s(N);
    for (size_t ia=0;ia<N;ia++)
        for (size_t ib=ia;ib<N;ib++)
            s(ia,ib)=ab->radials[ia]->Integrate(t2C,*ab->radials[ib],ab->pols[ia],ab->pols[ib],cl)*ab->ns[ia]*ab->ns[ib];

    return s;
}

} //namespace