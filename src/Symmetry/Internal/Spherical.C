// File: Symmetry/Internal/Imp/Spherical.C all spherical symmetry implementations
module;
#include <iosfwd>
export module qchem.Symmetry.Internal.Spherical;
export import qchem.Symmetry.Angular;

export namespace Symmetryns::Internal::Spherical
{
class Yl_Sym : public virtual SphericalSym
{
public:
    Yl_Sym(size_t l) : itsL(l) {};

    virtual size_t SequenceIndex() const {return itsL;} //Used for op<
    virtual size_t GetDegeneracy() const {return 2*itsL+1;}
    virtual size_t Getl         () const {return itsL;}
    virtual ivec_t Getmls       () const {return {};}
    virtual std::ostream&  Write(std::ostream&) const;
protected:
    size_t itsL;
    static const size_t LMAX=4;
};

class Ylm_Sym
    : public virtual SphericalSym
    , private Yl_Sym
{
public:
    Ylm_Sym(size_t l, const ivec_t& mls);

    virtual size_t SequenceIndex() const; //Used for op<
    virtual size_t GetDegeneracy() const;
    using Yl_Sym::Getl;
    virtual ivec_t Getmls() const {return mls;}

    virtual std::ostream&  Write(std::ostream&) const;
   
protected:
    ivec_t mls;
};

class Omega_k_Sym : public virtual SphericalSpinorSym
{
public:
    Omega_k_Sym(int κ);

    virtual size_t SequenceIndex() const; //Used for op<
    virtual size_t GetDegeneracy() const;
    virtual int    Getκ  () const {return κ;}
    virtual rvec_t Getmjs() const {return {};}
 
    virtual std::ostream&  Write(std::ostream&) const;

    
protected:
    int κ;
    static const size_t LMax=4;
};

class Omega_kmj_Sym
    : public virtual SphericalSpinorSym
    , private Omega_k_Sym

{
public:
    Omega_kmj_Sym(int κ, const rvec_t& mjs);

    virtual size_t SequenceIndex() const; //Used for op<
    virtual size_t GetDegeneracy() const;
    using Omega_k_Sym::Getκ;
    virtual rvec_t Getmjs() const {return mjs;}

    virtual std::ostream&  Write(std::ostream&) const;
   

protected:
    rvec_t mjs;
};

} //namespace