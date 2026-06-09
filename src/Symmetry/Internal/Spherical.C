// File: Symmetry/Internal/Imp/Spherical.C all spherical symmetry implementations
module;
#include <iosfwd>
export module qchem.Symmetry.Internal.Spherical;
export import qchem.Symmetry.Spherical;

export namespace Symmetry::Internal::Spherical
{
class Yl : public virtual Spherical
{
public:
    Yl(size_t l) : itsL(l) {};

    virtual size_t SequenceIndex() const {return itsL;} //Used for op<
    virtual size_t GetDegeneracy() const {return 2*itsL+1;}
    virtual size_t Getl         () const {return itsL;}
    virtual ivec_t Getmls       () const {return {};}
    virtual std::ostream&  Write(std::ostream&) const;
protected:
    size_t itsL;
    static const size_t LMAX=4;
};

class Ylm
    : public virtual Spherical
    , private Yl
{
public:
    Ylm(size_t l, const ivec_t& mls);

    virtual size_t SequenceIndex() const; //Used for op<
    virtual size_t GetDegeneracy() const;
    using Yl::Getl;
    virtual ivec_t Getmls() const {return mls;}

    virtual std::ostream&  Write(std::ostream&) const;
   
protected:
    ivec_t mls;
};

class Ωκ : public virtual SphericalSpinor
{
public:
    Ωκ(int κ);

    virtual size_t SequenceIndex() const; //Used for op<
    virtual size_t GetDegeneracy() const;
    virtual int    Getκ  () const {return κ;}
    virtual rvec_t Getmjs() const {return {};}
 
    virtual std::ostream&  Write(std::ostream&) const;

    
protected:
    int κ;
    static const size_t LMax=4;
};

class Ωκmj
    : public virtual SphericalSpinor
    , private Ωκ

{
public:
    Ωκmj(int κ, const rvec_t& mjs);

    virtual size_t SequenceIndex() const; //Used for op<
    virtual size_t GetDegeneracy() const;
    using Ωκ::Getκ;
    virtual rvec_t Getmjs() const {return mjs;}

    virtual std::ostream&  Write(std::ostream&) const;
   

protected:
    rvec_t mjs;
};

} //namespace