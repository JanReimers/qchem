// File: BasisSet/Molecule/PG_Spherical/Imp/Orbital_IBS.C  Spherical-Gaussian 2-centre fit integrals.
//
// The orbital 1E / 3C (DFT) / 4C (HF) integral building lives entirely in the Molecule-generic,
// evaluator-templated mixins (qchem.BasisSet.Molecule.IBS), instantiated with the spherical NR_Evaluator.
// Only the EFit_IBS 2-centre fit integrals stay here -- they go through the evaluator's transform-summed
// Overlap/Repulsion2C kernels (the spherical analog of PG_Cart's raw-radial loops).
module;
#include <cassert>

module qchem.BasisSet.Molecule.PG_Spherical;
import qchem.Blaze;

namespace BasisSet::Molecule::PG_Spherical
{

rsmat_t MakeOverlap2C(const Sph::NR_Evaluator* ab)
{
    assert(ab);
    size_t N=ab->size();
    rsmat_t s(N);
    for (size_t ia=0;ia<N;ia++)
        for (size_t ib=ia;ib<N;ib++)
            s(ia,ib)=ab->Overlap(ia,ib);
    return s;
}

rsmat_t MakeRepulsion2C(const Sph::NR_Evaluator* ab)
{
    assert(ab);
    size_t N=ab->size();
    rsmat_t s(N);
    for (size_t ia=0;ia<N;ia++)
        for (size_t ib=ia;ib<N;ib++)
            s(ia,ib)=ab->Repulsion2C(ia,ib);
    return s;
}

} //namespace
