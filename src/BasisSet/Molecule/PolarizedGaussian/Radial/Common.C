// File: Radial/Common.C  Partial implementation for the radial part of a basis function.
module;
#include <iosfwd>
export module qchem.BasisSet.Molecule.PolarizedGaussian.Internal.Radial.Common;
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.MnD.Hermite1;
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.RadialFunction;
import Common.UniqueIDImp;
import qchem.Types;

export namespace PolarizedGaussian
{

//
//  All radial function must store a list pairwise charge distributions.
//  They also all have a center and maximum orbital angular momentum L.
//  This commonality is implementated here.
//
class RadialCommon
    : public virtual RadialFunction
    , private UniqueIDImp
{
public:
    RadialCommon(                         );
    RadialCommon(const rvec3_t& Center,int L);
    RadialCommon(const RadialCommon&);
    ~RadialCommon(                        );

    virtual       bool      operator==(const RadialFunction&) const;
    virtual const rvec3_t&    GetCenter (                     ) const
    {
        return itsCenter;
    }
    virtual       int       GetL      (                     ) const
    {
        return itsL;
    }
    virtual const Hermite1& GetH1     (                     ) const;
    
    using UniqueIDImp::GetID;
    virtual std::ostream&           Write(std::ostream&) const;

protected:
    virtual Hermite1* MakeH1() const=0;

    rvec3_t             itsCenter;
    int               itsL;
    mutable Hermite1* itsH1;
};

} //namespace PolarizedGaussian


