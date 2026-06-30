// File: Symmetry/Internal/Imp/Spherical.C all spherical symmetry implementations
module;
#include <iosfwd>
export module qchem.Symmetry.Internal.Spherical;
export import qchem.Symmetry.Spherical;

export namespace qchem::Symmetry::Internal::Spherical
{
inline constexpr size_t LMax = 4;

class Yl : public virtual Spherical
{
public:
    Yl(size_t l) : itsL(l) {};

    virtual size_t SequenceIndex() const override {return itsL;} //Used for op<
    virtual size_t GetDegeneracy() const override {return 2*itsL+1;}
    virtual size_t Getl         () const override {return itsL;}
    virtual ivec_t Getmls       () const override {return {};}
    virtual std::ostream&  Write(std::ostream&) const override;
protected:
    size_t itsL;
};

class Ylm
    : public virtual Spherical
{
public:
    Ylm(size_t l, const ivec_t& mls);

    virtual size_t SequenceIndex() const override; //Used for op<
    virtual size_t GetDegeneracy() const override;
    virtual size_t Getl() const override {return itsL;}
    virtual ivec_t Getmls() const override {return mls;}

    virtual std::ostream&  Write(std::ostream&) const override;
   
protected:
    size_t itsL;
    ivec_t mls;
};

class Ωκ : public virtual SphericalSpinor
{
public:
    Ωκ(int κ);

    virtual size_t SequenceIndex() const override; //Used for op<
    virtual size_t GetDegeneracy() const override;
    virtual int    Getκ  () const override {return κ;}
    virtual rvec_t Getmjs() const override {return {};}
 
    virtual std::ostream&  Write(std::ostream&) const override;

    
protected:
    int κ;
};

class Ωκmj
    : public virtual SphericalSpinor
{
public:
    Ωκmj(int κ, const rvec_t& mjs);

    virtual size_t SequenceIndex() const override; //Used for op<
    virtual size_t GetDegeneracy() const override;
    virtual int    Getκ() const override {return κ;}
    virtual rvec_t Getmjs() const override {return mjs;}

    virtual std::ostream&  Write(std::ostream&) const override;
   

protected:
    int κ;
    rvec_t mjs;
};

} //namespace