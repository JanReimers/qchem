// File: BasisSet1/Atom/Evaluators/BS_Evaluator.C Interface for atom basis set evaluators.
module;
export module qchem.BasisSet1.Atom.Evaluators.BS;
export import qchem.BasisSet1.Atom.Evaluators.IBS;
export import qchem.BasisSet1.Internal.Cache4;
export import qchem.BasisSet1.Internal.ERI4;
import qchem.BasisSet1.Atom.Evaluators.Internal.AngularIntegrals;

export class BS_Evaluator : public Cache4
{
public:
    typedef AngularIntegrals::rvec11_t rvec11_t;
    virtual ~BS_Evaluator() {}; //g++ 15.2 BUG Compiler implemented destructor not created.

    virtual void Register(IBS_Evaluator*)=0; //Used for grouping all unique exponents or splines.
    virtual double loop_4_direct  (size_t id, size_t la, size_t lc,const rvec11_t& Ak) const=0; //Return vector dot product A[k]*R[k] 
    virtual double loop_4_exchange(size_t id, size_t la, size_t lc,const rvec11_t& Ak) const=0; //Return vector dot product A[k]*R[k] 
    // Angular integrals use Wigner-3j symbols and are basis set independent.
    rvec11_t Coulomb_AngularIntegrals(const IBS_Evaluator* a,const IBS_Evaluator* c) const;
    rvec11_t ExchangeAngularIntegrals(const IBS_Evaluator* a,const IBS_Evaluator* c) const;
    // These are the 4-nested hot loops, implemented at the general Atom level, using the virtual calls above.
    ERI4 Direct  (const IBS_Evaluator* a, const IBS_Evaluator* c) const;
    ERI4 Exchange(const IBS_Evaluator* a, const IBS_Evaluator* c) const;
};


