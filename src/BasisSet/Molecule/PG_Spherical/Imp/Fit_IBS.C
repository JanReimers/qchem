// File: BasisSet/Molecule/PG_Spherical/Imp/Fit_IBS.C  Spherical-Gaussian fit basis set, for MO calcs.
module;
#include <cassert>
#include <vector>

module qchem.BasisSet.Molecule.PG_Spherical;
import qchem.Blaze;

namespace BasisSet::Molecule::PG_Spherical
{

// Per-fit-component charge.  EFit_IBS IS-A spherical NR_Evaluator, so Charge(i) is the inherited
// transform-summed monopole (zero for the l>0 harmonics, nonzero only for the s components).
rvec_t EFit_IBS::MakeCharge() const
{
    rvec_t c(size());
    for (size_t i=0;i<size();i++) c[i]=Sph::NR_Evaluator::Charge(i);  // (Fit_IBS also declares Charge)
    return c;
}

// Cross fit-fit Coulomb repulsion (this fit x another fit), through the evaluator's cross Repulsion2C.
rmat_t EFit_IBS::MakeRepulsion(const FIT_CD_ABS& _b) const
{
    const Sph::NR_Evaluator* b=dynamic_cast<const Sph::NR_Evaluator*>(&_b);
    assert(b);
    size_t Na=size(),Nb=b->size();
    rmat_t s(Na,Nb);
    for (size_t ia=0;ia<Na;ia++)
        for (size_t ib=0;ib<Nb;ib++)
            s(ia,ib)=Repulsion2C(ia,*b,ib);
    assert(!blazem::isnan(s));
    return s;
}

} //namespace
