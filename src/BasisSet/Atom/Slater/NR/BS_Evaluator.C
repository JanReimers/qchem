// File: BasisSet/Atom/Slater/NR/BS_Evaluator.C
module;
export module BasisSet.Atom.Slater.NR.BS_Evaluator;
export import qchem.BasisSet.Atom.Gaussian.NR.BS_Evaluator;
import qchem.BasisSet.Atom.Internal.ExponentGrouper;
export import qchem.BasisSet.Atom.Rk;


export class Slater_BS 
    : public /* virtual g++-15.2 BUG failed to read compiled module cluster 32: Bad file data */ BS_Evaluator
{
public:
    virtual void Register(IBS_Evaluator *);
    virtual Rk*  Create (size_t ia,size_t ic,size_t ib,size_t id) const; //4 center
    virtual double loop_4_direct  (size_t id, size_t la, size_t lc,rvec11_t Ak) const; //Return vector dot product A[k]*R[k] 
    virtual double loop_4_exchange(size_t id, size_t la, size_t lc,rvec11_t Ak) const; //Return vector dot product A[k]*R[k] 
private:
    ExponentGrouper grouper;
};

