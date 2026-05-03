// File: IBS_Common.C  Irrep Basis set common implementation.
module;
#include <iosfwd>
#include "forward.H"
export module qchem.BasisSet1.Internal.Common;
export import qchem.Types;
import Common.UniqueIDImp;
export import qchem.BasisSet1;


export class BS_Common1
    : public virtual BasisSet1
    , private UniqueIDImp
{
public:
    BS_Common1() {};
    virtual ~BS_Common1() {};

    virtual size_t GetNumFunctions() const;
    virtual irrepv_t GetIrreps(const Spin& ms) const;
//
//  StreamableObject stuff.
//
    virtual std::ostream&  Write(std::ostream&    ) const;
protected:
    void Insert(Orbital_1E_IBS1<double>*);

    virtual const_iterator begin() const {return itsBasisSets.begin();}
    virtual const_iterator end  () const {return itsBasisSets.end  ();}
    friend class SlaterRadialIntegralTests;
    friend class DiracIntegralTests;

    bsv_t itsBasisSets;
};

