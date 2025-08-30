// File: BasisSet/Atom/radial/Gaussian_BS.C
module;
export module BasisSet.Atom.Gaussian_BS;
export import BasisSet.Atom.BS_Evaluator;
export import qchem.BasisSet.Atom.Rk;
import qchem.BasisSet.Atom.Internal.ExponentGrouper;


export class Gaussian_BS 
    : public virtual BS_Evaluator
{
public:
    virtual void Register(IBS_Evaluator *);
    virtual Rk*  Create (size_t ia,size_t ic,size_t ib,size_t id) const; //4 center
    virtual RVec loop_4_direct  (size_t id, size_t la, size_t lc) const;
    virtual RVec loop_4_exchange(size_t id, size_t la, size_t lc) const;
private:
    ExponentGrouper grouper;
};

