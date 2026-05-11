// File: BasisSet1/Atom/Evaluators/Gaussian/BS_Evaluator.C
module;
export module qchem.BasisSet1.Atom.Evaluators.Gaussian.BS; 
export import qchem.BasisSet1.Atom.Evaluators.Gaussian.IBS; 
export import qchem.BasisSet1.Atom.Evaluators.BS;
export import qchem.BasisSet1.Atom.Rk;
import qchem.BasisSet1.Atom.Internal.ExponentGrouper;


export class Gaussian_BS_Evaluator 
    : public /* virtual g++-15.2 BUG failed to read compiled module cluster 32: Bad file data */ BS_Evaluator
{
public:
    using IBS_Evaluator_t = Gaussian_IBS_Evaluator;
    virtual void Register(IBS_Evaluator *);
    void BuildCache(size_t LMax) {};
    virtual Rk*  Create (size_t ia,size_t ic,size_t ib,size_t id) const; //4 center
    virtual double loop_4_direct  (size_t id, size_t la, size_t lc,const rvec11_t& Ak) const; //Return vector dot product A[k]*R[k] 
    virtual double loop_4_exchange(size_t id, size_t la, size_t lc,const rvec11_t& Ak) const; //Return vector dot product A[k]*R[k] 

private:
    ExponentGrouper grouper;
};

