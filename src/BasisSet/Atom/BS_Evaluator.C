// File: BasisSet/Atom/BS_Evaluator.C Create Rk structures for 2 electron repulsion integrals (2ERIs).
module;

export module BasisSet.Atom.BS_Evaluator;
export import qchem.BasisSet.Atom.Internal.ExponentGrouper;
export import qchem.BasisSet.Atom.internal.Rk;

export class BS_Evaluator
{
public:
    virtual void Register(ExponentGrouper*)=0; //Set up unique spline or exponent indexes.
     virtual Rk* CreateRk (size_t ia,size_t ic,size_t ib,size_t id) const=0; //4 center
};