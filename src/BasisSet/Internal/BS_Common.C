// File: IBS_Common.C  Irrep Basis set common implementation.
module;
#include <iosfwd>
class SlaterRadialIntegralTests;
class DiracIntegralTests;

export module qchem.BasisSet.Internal.Common;
export import qchem.Types;
import qchem.LAParams;
import Common.UniqueIDImp;
import qchem.BasisSet;


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
    void Insert(Orbital_IBS<double>*);

    virtual const_iterator begin() const {return itsBasisSets.begin();}
    virtual const_iterator end  () const {return itsBasisSets.end  ();}
    friend class SlaterRadialIntegralTests;
    friend class DiracIntegralTests;

    bsv_t itsBasisSets;
};

