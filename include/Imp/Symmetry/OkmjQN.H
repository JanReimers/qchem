// File: Symmetry/OkmjQN.H  Spherical Spinor Omega_kmj symmetry.
#ifndef _Symmetry_OkmjQN_H_
#define _Symmetry_OkmjQN_H_

#include <QuantumNumber.H>
#include <utility>


//---------------------------------------------------------------------------------
//
//  Spherical Spinor with total AM J=L+S.  QNs are kappa and mj.
//
class Omega_kmjQN
    : public virtual QuantumNumber

{
public:
    Omega_kmjQN(                 );
    Omega_kmjQN(int kappa, double mj);

    virtual bool Match(const QuantumNumber&) const;
    virtual int GetDegeneracy() const;
    QuantumNumber* AddPrincipleQN(int index) const;
    virtual std::pair<int,int> GetN(const int (&N)[4], const int (&Nv)[4], int NUnpaired) const;

    virtual std::ostream&       Write(std::ostream&) const;
    virtual std::istream&       Read (std::istream&)      ;
    virtual QuantumNumber* Clone(        ) const;

    int    GetKappa() const {return kappa;}
    double Getj    () const {return j(kappa);}
    int    Getl    () const {return l(kappa);}
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
