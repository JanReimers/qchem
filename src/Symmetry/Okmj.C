// File: Symmetry/Okmj.C  Spherical Spinor Omega_kmj symmetry.
module;
#include <iosfwd>
#include <utility> //std::pair
#include <vector>

export module qchem.Symmetry.Okmj;
export import qchem.Symmetry.Angular;
export import qchem.Types;
//---------------------------------------------------------------------------------
//
//  Spherical Spinor with total AM J=L+S.  QNs are kappa (mj gets averaged over).
//
export class Omega_k_Sym
    : public virtual Angular_Sym

{
public:
    Omega_k_Sym(int kappa);

    virtual size_t SequenceIndex() const; //Used for op<
    virtual int GetDegeneracy() const;
    
    virtual std::ostream&  Write(std::ostream&) const;

    int    GetKappa() const {return kappa;}
    double Getj    () const {return j(kappa);}
    int    GetL    () const {return l(kappa);}

    static double j(int kappa) {return kappa>0 ? kappa-0.5 : -kappa-0.5;}
    static double l(int kappa) {return kappa>0 ? kappa : -kappa-1;}
protected:
    int kappa;
    static const int LMax=4;
};


//---------------------------------------------------------------------------------
//
//  Spherical Spinor with total AM J=L+S.  QNs are kappa and mj.
//
export class Omega_kmj_Sym
    : public virtual Angular_Sym

{
public:
    Omega_kmj_Sym(int kappa, const std::vector<double>& mjs);

    virtual size_t SequenceIndex() const; //Used for op<
    virtual int GetDegeneracy() const;
    Symmetry* AddPrincipleQN(int index) const;

    virtual std::ostream&  Write(std::ostream&) const;
   
    int    GetKappa() const {return kappa;}
    double Getj    () const {return j(kappa);}
    int    GetL    () const {return l(kappa);}

    static double j(int kappa) {return kappa>0 ? kappa-0.5 : -kappa-0.5;}
    static double l(int kappa) {return kappa>0 ? kappa : -kappa-1;}
    static double ms(int kappa) {return kappa<0 ? -0.5 : 0.5;}
    static int    ml(int kappa, double mj) {return mj-ms(kappa);}
protected:
    int kappa;
    std::vector<double> mjs;
    static const int LMax=3;
};

