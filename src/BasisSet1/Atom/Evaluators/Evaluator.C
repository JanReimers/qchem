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
    Evaluator(int _l) : l(_l), grouper(0) {};
    Evaluator(const sym_t& ylm);
    virtual ~Evaluator() {};

    virtual void          Register     (Grouper*)=0; //Set up unique spline or exponent indexes.
    virtual size_t        size         () const = 0;
    virtual size_t        maxSpan      () const {return size();}  //assume no overlap for indeces separated by > maxSpan
    virtual size_t        GetVectorSize() const {return size();}
    virtual int           Getl         () const {return l;}
    virtual       rvec_t  Norm         () const = 0;

    iv_t                  indices      (             ) const {return iv_t(size_t(0),size());}
    iv_t                  indices      (size_t start ) const {return iv_t(start,size());}
    virtual size_t        es_index     (size_t i     ) const {return es_indices[i];}

    virtual rvec11_t      CoulombAk (const Evaluator& other) const = 0;
    virtual rvec11_t      ExchangeAk(const Evaluator& other) const = 0;

    virtual std::ostream& Write  (std::ostream&) const=0;
    virtual std::string RadialID () const=0;
    virtual std::string AngularID() const;
    virtual std::string Name     () const=0;

protected:
    friend class Cache4Tests;

    int    l;
    const  ExponentGrouper* grouper;
    std::vector<size_t> es_indices; //Unique exponent index

};




} //namespace