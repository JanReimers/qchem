// File: Symmetry/OkmjQN.H  Spherical Spinor Omega_kmj symmetry.
#ifndef _Symmetry_OkmjQN_H_
#define _Symmetry_OkmjQN_H_

#include "Imp/Symmetry/AngularQN.H"
#include <utility>

//---------------------------------------------------------------------------------
//
//  Spherical Spinor with total AM J=L+S.  QNs are kappa (mj gets averaged over).
//
class Omega_kQN
    : public virtual AngularQN

{
public:
    Omega_kQN(         );
    Omega_kQN(int kappa);

    virtual bool Match(const QuantumNumber&) const;
    virtual int GetDegeneracy() const;
    QuantumNumber* AddPrincipleQN(int index) const;
    virtual std::pair<int,int> GetN(const int (&N)[4], const int (&Nv)[4], int NUnpaired) const;

    virtual std::ostream&  Write(std::ostream&) const;
    virtual std::istream&  Read (std::istream&)      ;
    virtual AngularQN* Clone(        ) const;

    int    GetKappa() const {return kappa;}
    double Getj    () const {return j(kappa);}
    int    GetL    () const {return l(kappa);}

    static double j(int kappa) {return kappa>0 ? kappa-0.5 : -kappa-0.5;}
    static double l(int kappa) {return kappa>0 ? kappa : -kappa-1;}
protected:
    int kappa;
    static const int LMax=3;
};


//---------------------------------------------------------------------------------
//
//  Spherical Spinor with total AM J=L+S.  QNs are kappa and mj.
//
class Omega_kmjQN
    : public virtual AngularQN

{
public:
    Omega_kmjQN(                 );
    Omega_kmjQN(int kappa, double mj);

    virtual bool Match(const QuantumNumber&) const;
    virtual int GetDegeneracy() const;
    QuantumNumber* AddPrincipleQN(int index) const;
    virtual std::pair<int,int> GetN(const int (&N)[4], const int (&Nv)[4], int NUnpaired) const;

    virtual std::ostream&  Write(std::ostream&) const;
    virtual std::istream&  Read (std::istream&)      ;
    virtual AngularQN* Clone(        ) const;

    int    GetKappa() const {return kappa;}
    double Getj    () const {return j(kappa);}
    int    GetL    () const {return l(kappa);}
    double Getmj   () const {return mj;}
    int    Getml   () const {return ml(kappa,mj);}

    static double j(int kappa) {return kappa>0 ? kappa-0.5 : -kappa-0.5;}
    static double l(int kappa) {return kappa>0 ? kappa : -kappa-1;}
    static double ms(int kappa) {return kappa<0 ? -0.5 : 0.5;}
    static int    ml(int kappa, double mj) {return mj-ms(kappa);}
protected:
    virtual std::pair<int,int> GetNk(const int (&N)[4], const int (&Nv)[4], int NUnpaired) const;

    int kappa;
    double mj;
    static const int LMax=3;
};


#endif //_Symmetry_OkmjQN_H_
