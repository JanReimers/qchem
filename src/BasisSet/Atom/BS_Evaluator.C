// File: BasisSet/Atom/BS_Evaluator.C Generic hot loop 2 electron repulsion integrals (2ERIs).
module;
export module qchem.BasisSet.Atom.BS_Evaluator;
export import qchem.BasisSet.Atom.IBS_Evaluator;
export import qchem.BasisSet.Internal.Cache4;
export import qchem.BasisSet.Internal.ERI4;


export class BS_Evaluator : public Cache4
{
public:
    virtual void Register(IBS_Evaluator*)=0; //Used for grouping all unique exponents or splines.
    virtual RVec loop_4_direct  (size_t id, size_t la, size_t lc) const=0; //Return a vector of R[k] values 
    virtual RVec loop_4_exchange(size_t id, size_t la, size_t lc) const=0; //Return a vector of R[k] values 
    // Angular integrals use Wigner-3j symbols and are basis set independent.
    RVec Coulomb_AngularIntegrals(const IBS_Evaluator* a,const IBS_Evaluator* c) const;
    RVec ExchangeAngularIntegrals(const IBS_Evaluator* a,const IBS_Evaluator* c) const;
    // These are the 4-nested hot loops, implemented at the general Atom level, using the virtual calls above.
    ERI4 Direct  (const IBS_Evaluator* a, const IBS_Evaluator* c) const;
    ERI4 Exchange(const IBS_Evaluator* a, const IBS_Evaluator* c) const;
};


