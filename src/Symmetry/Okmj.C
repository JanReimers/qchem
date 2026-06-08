// File: Symmetry/Okmj.C  Spherical Spinor Omega_kmj symmetry.
module;
#include <iosfwd>
#include <utility> //std::pair
#include <vector>

export module qchem.Symmetry.Okmj;
export import qchem.Symmetry.Angular;
// export import qchem.Types;
//---------------------------------------------------------------------------------
//
//  Spherical Spinor with total AM J=L+S.  QNs are kappa (mj gets averaged over).
//
export class Omega_k_Sym : public virtual SphericalSpinorSym
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


//---------------------------------------------------------------------------------
//
//  Spherical Spinor with total AM J=L+S.  QNs are kappa and mj.
//
export class Omega_kmj_Sym
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

