// File: BasisSet/Molecule/PG_Cart/Imp/Fit_IBS.C  Polarized Gaussian fit basis set, for MO calculations.
module;
#include <cassert>
// #include <algorithm> //Need std::max
#include <iostream>
#include <fstream>
#include <string>
#include <memory>
#include <vector>

module qchem.BasisSet.Molecule.PG_Cart;
import qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.GaussianRF;
import qchem.BasisSet.Molecule.Readers.Gaussian94;
import qchem.BasisSet;
import qchem.Cluster;
import qchem.Symmetry.Unit;
import qchem.stl_io;
import qchem.Streamable;
import qchem.Math;
import qchem.Blaze;

namespace BasisSet::Molecule::PG_Cart
{
using namespace ::BasisSet::Molecule::Evaluators::PG_Cart_MnD;  // Cartesian glue moved out to PG_Cart_MnD

rvec_t EFit_IBS::MakeCharge() const
{
    const PGData* a=dynamic_cast<const PGData*>(this);
    assert(a);
    rvec_t c(a->size());
    int i=0;
    for (auto r:a->radials)
    {
        c[i]=r->GetCharge(a->pols[i])*a->ns[i]; 
        i++;       
    }

    // 1 line copilot version
    //for (auto i:ab->ns.indices())  c(i)=ab->radials[i-1]->GetCharge(ab->pols[i-1])*ab->ns(i);
    return c;
}

rmat_t EFit_IBS::MakeRepulsion(const Fit_IBS& _b) const
{   
    const PGData* a=dynamic_cast<const PGData*>(this);
    const PGData* b=dynamic_cast<const PGData*>(&_b);
    assert(a);
    assert(b); 
    int Na=a->size(),Nb=b->size();
    rmat_t s(Na,Nb);
    for (size_t ia=0;ia<Na;ia++)
        for (size_t ib=0;ib<Nb;ib++)
            s(ia,ib)=a->radials[ia]->Repulsion2C(
                *b->radials[ib],a->pols[ia],b->pols[ib])*a->ns[ia]*b->ns[ib];
    assert(!blazem::isnan(s));
    return s;
}

} //namespace