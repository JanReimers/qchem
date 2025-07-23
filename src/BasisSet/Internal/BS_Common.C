// File: IBS_Common.C  Irrep Basis set common implementation.
module;
#include <iosfwd>
class SlaterRadialIntegralTests;
class DiracIntegralTests;

export module qchem.BasisSet.Internal.Common;
import qchem.LAParams;
import Common.UniqueIDImp;
import qchem.BasisSet;
export import qchem.Types;


export class BS_Common
    : public virtual BasisSet
    , private UniqueIDImp
{
public:
    BS_Common() {};
    virtual ~BS_Common() {};
    virtual void Set(const LAParams&);

    virtual size_t GetNumFunctions() const;
    virtual symv_t GetSymmetries  () const;

//
//  StreamableObject stuff.
//
    virtual std::ostream&  Write(std::ostream&    ) const;
protected:
    typedef TOrbital_IBS<double> bs_t; 
    virtual void Insert(bs_t*);

    virtual const_iterator begin() const {return itsBasisSets.begin();}
    virtual const_iterator end  () const {return itsBasisSets.end  ();}
    friend class SlaterRadialIntegralTests;
    friend class DiracIntegralTests;

    bsv_t itsBasisSets;
};

