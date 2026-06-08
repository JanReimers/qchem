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

import qchem.BasisSet.Internal.Cache4;

export namespace BasisSet::Atom::Evaluators
{
using dERI3=ERI3<double>;


class Evaluator 
    : public virtual Cache4_Client
    , public VectorFunction<double>
{
public:
    Evaluator(int _l, const ivec_t& _mls) : l(_l), mls(_mls),ns(0),grouper(0) {};
    Evaluator(const sym_t& ylm);
    virtual ~Evaluator() {};

    virtual void          Register     (Grouper*)=0; //Set up unique spline or exponent indexes.
    virtual size_t        size         () const {return ns.size();}
    virtual size_t        maxSpan      () const {return size();}  //assume no overlap for indeces separated by > maxSpan
    virtual size_t        GetVectorSize() const {return size();}
    virtual int           Getl         () const {return l;}
    virtual const ivec_t& Getmls       () const {return mls;}
    virtual       rvec_t  Norm         () const {return ns;}

    iv_t                  indices      (             ) const {return iv_t(size_t(0),size());}
    iv_t                  indices      (size_t start ) const {return iv_t(start,size());}
    virtual size_t        es_index     (size_t i     ) const {return es_indices[i];}

   
    virtual std::ostream& Write  (std::ostream&) const=0;
    virtual std::string RadialID () const=0;
    virtual std::string AngularID() const;
    virtual std::string Name     () const=0;

    static rvec11_t Coulomb_AngularIntegrals(const Evaluator& a,const Evaluator& c);
    static rvec11_t ExchangeAngularIntegrals(const Evaluator& a,const Evaluator& b);
protected:
    friend class Cache4Tests;

    int    l;
    ivec_t mls;
    rvec_t ns;
    const  ExponentGrouper* grouper;
    std::vector<size_t> es_indices; //Unique exponent index

};




} //namespace