// File: RadialCommon.H  Partial implementation for the radial part of a basis function.
#ifndef _RadialCommon_H_
#define _RadialCommon_H_



#include "../RadialFunction.H"
#include "Imp/Misc/UniqueID/UniqueIDImp.H"

namespace PolarizedGaussian
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
    RadialCommon(const RVec3& Center,int L);
    RadialCommon(const RadialCommon&);
    ~RadialCommon(                        );

    virtual       bool      operator==(const RadialFunction&) const;
    virtual const RVec3&    GetCenter (                     ) const
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
    virtual std::istream&           Read (std::istream&)      ;

private:
    virtual Hermite1* MakeH1() const=0;

    RVec3             itsCenter;
    int               itsL;
    mutable Hermite1* itsH1;
};

} //namespace PolarizedGaussian

#endif //_RadialCommon_H_
