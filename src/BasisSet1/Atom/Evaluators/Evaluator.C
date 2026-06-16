// File: BasisSet/Atom/Evaluator.C
module;
#include <vector>
#include <iosfwd>
#include "forward.H"
export module qchem.BasisSet.Atom.Evaluators;
export import qchem.BasisSet.Atom.Evaluators.Internal.ExponentGrouper;
export import qchem.VectorFunction;
export import qchem.Streamable;

import qchem.BasisSet.Internal.Cache4;

export namespace BasisSet::Atom::Evaluators
{
//
//  Abstract inteface for all evaluators.  We need to use mixins when creating the final (concrete) evaluators in order to support
//  Any combination of:
//      {Radial basis: BS,SG,SL}*{Physics: NR(radial)+Spherical(Angluar),RKB(radial)+SphericalSpiner(angular)}*{Hamiltonian: 1E,DFT,HF}
//  This gives a total of 3x2x3=18 combinations.  In practice it is not that bad. For each of the 6 radial/physics combinations we
//  just need final evaluators that support as much as possible out of set {1E,DFT,HF}.
//  The mixin flexiblity is important when working a new basis set type.  You ar not force the implement everything required for all
//  18 combinations.  A good starting point is the just implement Physics: NR(radial)+Spherical(Angluar) for Hamiltonian:1E and then add
//  support for more integrals later. Spherical(Angluar) is independent Radial basis so it is already implemented.
//
//  We use these abstract interfaces to define what virtual functions are required for various evaluators.  These virtual
//    functions only get used outside hot loops.
//  We use concepts (./Concepts.C) to define what inline (non-virtual) functions are required.  These mostly used
//    in hot loops.
//
//  This mixture of OOD methodolgies obviously adds complexity, but is required for efficiency and good design.
//
class Evaluator
    : public virtual Streamable
    , public VectorFunction<double> //Try virtual and get: virtual function 'VectorFunction<double>::GetVectorSize' has more than one final overrider in ...
{
public:
    virtual size_t  GetVectorSize() const {return size();}
    
    virtual int     Getl    () const = 0; //  l is Used **everywhere**
    virtual size_t  size    () const = 0;
    virtual size_t  maxSpan () const {return size();}  //assume no overlap for indices separated by > maxSpan
    virtual rvec_t  Norm    () const = 0; //Normalization constants.
    virtual std::ostream& Write    (std::ostream&) const=0;
    virtual std::string   RadialID () const=0; // For creating keys into cache database.
    virtual std::string   Name     () const=0;
    // Helper functions for range based loops.
    iv_t indices(             ) const {return iv_t(size_t(0),size());} 
    iv_t indices(size_t start ) const {return iv_t(start    ,size());} //Use this for 2nd index of symmetric matrices.
};

class HF_Evaluator : public virtual Cache4_Client
{
public:
    virtual void    Register(Grouper*)=0;
    virtual Cache4* MakeCache4() const=0;
};

class Angular
{
public:
    virtual std::string AngularID () const=0; // For creating keys into cache database.
    // These are only required for HF and DHF ERI calculations. And fully implemented for HF and DHF.
    virtual rvec11_t DirectAk  (const Evaluator& other) const = 0;
    virtual rvec11_t ExchangeAk(const Evaluator& other) const = 0;
};

} //namespace