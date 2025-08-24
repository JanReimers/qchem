// File: BasisSet/Atom/BS_Evaluator.C Create Rk structures for 2 electron repulsion integrals (2ERIs).
module;
export module BasisSet.Atom.BS_Evaluator;
export import BasisSet.Atom.IBS_Evaluator;
export import qchem.BasisSet.Internal.Cache4;
import qchem.BasisSet.Atom.Internal.AngularIntegrals;

export class BS_Evaluator : public Cache4
{
public:
    virtual ~BS_Evaluator() {};
    virtual void Register(IBS_Evaluator*)=0;
    virtual RVec Coulomb_AngularIntegrals(const IBS_Evaluator* a,const IBS_Evaluator* c) const;
    virtual RVec ExchangeAngularIntegrals(const IBS_Evaluator* a,const IBS_Evaluator* c) const;
    virtual RVec loop_4_direct  (size_t id, size_t la, size_t lc) const=0;
    virtual RVec loop_4_exchange(size_t id, size_t la, size_t lc) const=0;
};

RVec BS_Evaluator::Coulomb_AngularIntegrals(const IBS_Evaluator* a,const IBS_Evaluator* c) const
{
    return AngularIntegrals::Coulomb(a->Getl(),c->Getl());
}

RVec BS_Evaluator::ExchangeAngularIntegrals(const IBS_Evaluator* a,const IBS_Evaluator* c) const
{
    return AngularIntegrals::Exchange(a->Getl(),c->Getl());
}