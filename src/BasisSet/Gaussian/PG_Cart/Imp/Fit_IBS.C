// File: BasisSet/Gaussian/PG_Cart/Imp/Fit_IBS.C  Polarized Gaussian fit basis set, for MO calculations.
module;
#include <cassert>
// #include <algorithm> //Need std::max
#include <iostream>
#include <fstream>
#include <string>
#include <memory>
#include <vector>
#include <stdexcept>

module qchem.BasisSet.Gaussian.PG_Cart;
import qchem.Blaze;   // blazem:: in the asserts: invisible to an NDEBUG build
import qchem.BasisSet.Gaussian.Evaluators.PG_Cart_MnD.GaussianRF;
import qchem.BasisSet.Gaussian.Evaluators.PG_LibCint;   // the Coulomb metric's derivative integrals (int2c2e_ip1), gated in tests/M_LibCint.C
import qchem.BasisSet;
import qchem.Structure;
import qchem.Streamable;

namespace qchem::BasisSet::Gaussian::PG_Cart
{
using namespace ::qchem::BasisSet::Gaussian::Evaluators::PG_Cart_MnD;  // Cartesian glue moved out to PG_Cart_MnD

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

rmat_t EFit_IBS::MakeRepulsion(const rFIT_CD_ABS& _b) const
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

Outcome<rvec3vec_t,AD::NoGradient> EFit_IBS::MakeRepulsionPullback(const rsmat_t& A, const Structure* st) const
{
    if (!st) throw std::invalid_argument("PG_Cart::EFit_IBS: a geometry pullback needs the structure");
    const PGData* a=dynamic_cast<const PGData*>(this);   // the block IS its PGData (as MakeCharge reads it)
    if (!a) throw std::logic_error("PG_Cart::EFit_IBS: the fit block carries no PGData");
    if (A.rows()!=a->size())
        throw std::invalid_argument("PG_Cart::EFit_IBS::MakeRepulsionPullback: the adjoint is "+std::to_string(A.rows())+" square, the fit has "
                                    +std::to_string(a->size())+" functions");
    const std::vector<rvec3_t> g=Evaluators::PG_LibCint::NR_Evaluator(*a, st).RepulsionTwoC_Gradient(A);
    rvec3vec_t F(g.size()); for (size_t B=0;B<g.size();B++) F[B]=g[B];
    return Outcome<rvec3vec_t,AD::NoGradient>::Ok(F);
}

} //namespace