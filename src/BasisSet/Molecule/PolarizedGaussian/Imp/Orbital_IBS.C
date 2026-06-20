// File: BasisSet/Molecule/PolarizedGaussian/Imp/Orbital_IBS.C  Polarized Gaussian fit basis set, for MO calculations.
module;
#include <cassert>
#include <algorithm> //Need std::max
#include <iostream>
#include <fstream>
#include <string>
#include <memory>
#include <vector>

module qchem.BasisSet.Molecule.PolarizedGaussian;
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.GaussianRF;
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.Readers.Gaussian94;
import qchem.BasisSet.Molecule.Evaluators;       // generic 1E matrix builders
import qchem.BasisSet.Molecule.Evaluators.PG;    // PG_Evaluator
import qchem.BasisSet.Internal.IntegralEnums;
import qchem.BasisSet;
import qchem.Cluster;
import qchem.Symmetry.Unit;
import qchem.stl_io;
import qchem.Streamable;
import qchem.Math;
import qchem.Blaze;


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

// The orbital 1E integrals now flow through the molecular Evaluator: wrap this IBS (an Orbital_IBS
// IS-A PGData) in a PG_Evaluator view and let the generic, basis-agnostic builders run the i,j loop.
// (EFit_IBS still uses MakeIntegrals for its Overlap2C/Repulsion2C fit integrals -- not 1E concepts.)
// MakeKinetic returns the <p^2>=<-nabla^2> block (no 1/2; see BasisSet/Orbital_1E_IBS.C).
rsmat_t Orbital_IBS::MakeOverlap() const
{
    return Evaluators::OverlapMatrix(Evaluators::PG_Evaluator(*this, nullptr));
}
rsmat_t Orbital_IBS::MakeKinetic() const
{
    return Evaluators::KineticMatrix(Evaluators::PG_Evaluator(*this, nullptr));
}
rsmat_t Orbital_IBS::MakeNuclear(const Cluster* cl) const
{
    return Evaluators::NuclearMatrix(Evaluators::PG_Evaluator(*this, cl));
}

// 3-centre (DFT) integrals through the evaluator: for each fit component ic, build the symmetric
// (ia,ib) block via the evaluator's ThreeC kernel (which folds in all three normalizations).
static ERI3<double> Make3C(qchem::IType3C type, const PGData& a, const PGData& c)
{
    Evaluators::PG_Evaluator aE(a, nullptr), cE(c, nullptr);
    size_t Na=aE.size(), Nc=cE.size();
    ERI3<double> s3;
    for (size_t ic=0; ic<Nc; ic++)
    {
        rsmat_t s(Na);
        for (size_t ia=0; ia<Na; ia++)
            for (size_t ib=ia; ib<Na; ib++)
                s(ia,ib)=aE.ThreeC(type, ia, aE, ib, cE, ic);
        s3.push_back(s);
    }
    return s3;
}
ERI3<double> Orbital_IBS::MakeOverlap3C(const Fit_IBS& _c) const
{
    return Make3C(qchem::Overlap3C, *this, *dynamic_cast<const PGData*>(&_c));
}
ERI3<double> Orbital_IBS::MakeRepulsion3C(const Fit_IBS& _c) const
{
    return Make3C(qchem::Repulsion3C, *this, *dynamic_cast<const PGData*>(&_c));
}

// 4-centre HF Coulomb (ab|cd): a,b on this orbital basis, c,d on the partner.  The intricate ERI loop
// + symmetry packing is unchanged (hoisting it is plan Goal D); only the per-element integral now goes
// through the evaluator's FourC kernel (which folds in all four normalizations).
ERI4 Orbital_IBS::MakeDirect  (const Orbital_HF_IBS<double>& _c) const
{
    const PGData* c=dynamic_cast<const PGData* >(&_c);
    assert(c);
    Evaluators::PG_Evaluator aE(*this, nullptr), cE(*c, nullptr);
    size_t Na=aE.size(), Nc=cE.size();
    ERI4 J(Na,Nc);

    for (size_t ia:iv_t(0,Na))
        for (size_t ib:iv_t(ia,Na))
        {
            rsmat_t& Jab=J(ia,ib);
            for (size_t ic:iv_t(0,Nc))
                for (size_t id:iv_t(ic,Nc))
                    Jab(ic,id)=aE.FourC(ia, aE, ib, cE, ic, cE, id);   // (a a | c c) slots
        }
    return J;
}

// 4-centre HF Exchange: slots (a b | a b).  Symmetry packing preserved exactly; element via FourC.
ERI4 Orbital_IBS::MakeExchange(const Orbital_HF_IBS<double>& _b) const
{
    const PGData* b=dynamic_cast<const PGData* >(&_b);
    assert(b);
    Evaluators::PG_Evaluator aE(*this, nullptr), bE(*b, nullptr);
    size_t Na=aE.size(), Nb=bE.size();
    ERI4 K(Na,Nb);
    for (size_t ia:iv_t(0,Na))
        for (size_t ib:iv_t(0,Nb))

            for (size_t ic:iv_t(ia,Na))
            {
                rsmat_t& Kac=K(ia,ic);
                for (size_t id:iv_t(0,Nb))
                {
                    double v=aE.FourC(ia, bE, ib, aE, ic, bE, id);   // (a b | a b) slots
                    if (ib==id)      Kac(ib,id) = v;
                    else if (ib<id)  Kac(ib,id)+= 0.5*v;
                    else             Kac(id,ib)+= 0.5*v;
                }
            }
    return K;
}


} //namespace