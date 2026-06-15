// File: BasisSet/Atom/Evaluator.C
module;
#include <vector>
#include <iosfwd>
#include "forward.H"
export module qchem.BasisSet.Atom.Evaluators.IBS;
export import qchem.BasisSet.Atom.Evaluators.Internal.ExponentGrouper;
export import qchem.BasisSet.Internal.ERI3;
export import qchem.Symmetry.Irrep;
export import qchem.VectorFunction;
export import qchem.Streamable;

import qchem.BasisSet.Internal.Cache4;

export namespace BasisSet::Atom::Evaluators
{
using dERI3=ERI3<double>;

//
//  Abstract inteface for all evaluators.  We need to use mixins when creating the final (concrete) evaluators in order to support
//  Any combination of {Radial basis: BS,SG,SL}*{Physics: NR(radial)+Spherical(Angluar),RKB(radial)+SphericalSpiner(angular)}*{1E,DFT,HF}
//  This gives a total of 3x2x3=18 combinations.  In practice it is not that bad. For each of the 6 radial/physics combinations we
//  just need final evaluators that support as much as possible out of set {1E,DFT,HF}.
//  The mixin flexiblity is important when working a new basis set type.  You ar not force the implement everything required for all
//  18 combinations.  A good starting point is the just implement 
//
class Evaluator
    : public virtual Cache4_Client
    , public virtual Streamable
    , public VectorFunction<double> //Try virtual and get: virtual function 'VectorFunction<double>::GetVectorSize' has more than one final overrider in ...
{
public:
    virtual size_t        GetVectorSize() const {return size();}
    //  Used **everywhere***
    virtual int     Getl    () const = 0;
    //  For radial
    virtual void    Register(Grouper*)=0; //Set up unique spline or exponent indexes.
    virtual size_t  size    () const = 0;
    virtual size_t  maxSpan () const {return size();}  //assume no overlap for indeces separated by > maxSpan
    virtual rvec_t  Norm    () const = 0;
            iv_t    indices (             ) const {return iv_t(size_t(0),size());} //For range loops
            iv_t    indices (size_t start ) const {return iv_t(start    ,size());} //For range loops
    virtual size_t  es_index(size_t i     ) const =0; //Get the index of basis function i in the grouper.
    // For angular

    // For Streamable.
    virtual std::ostream& Write  (std::ostream&) const=0;
    // For creating keys into cache database.
    virtual std::string RadialID () const=0;
    virtual std::string Name     () const=0;
};

class Angular
{
public:
    // For creating keys into cache database.
    virtual std::string AngularID () const=0;
    // These are only required for HF and DHF ERI calculations.
    virtual rvec11_t DirectAk  (const Evaluator& other) const = 0;
    virtual rvec11_t ExchangeAk(const Evaluator& other) const = 0;
};

} //namespace