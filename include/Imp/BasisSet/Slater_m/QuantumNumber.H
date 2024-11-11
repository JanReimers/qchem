// File: Slater_m/QuantumNumer.H  Spherical harmonic Ylm symmetry.
#ifndef _Slater_m_QuantumNumer_H_
#define _Slater_m_QuantumNumer_H_



#include "Imp/BasisSet/SphericalGaussian/QuantumNumber.H"


//---------------------------------------------------------------------------------
//
//  angular momentum L and m.
//

class YlmQN
    : public virtual QuantumNumber
    , public SphericalSymmetryQN
{
public:
    YlmQN(            );
    YlmQN(int l, int m);

    virtual bool Match(const QuantumNumber&) const;
    virtual int GetDegeneracy() const;
    virtual std::pair<int,int> GetN(const int (&N)[4], const int (&Nv)[4], int NUnpaired) const;

    virtual std::ostream&       Write(std::ostream&) const;
    virtual std::istream&       Read (std::istream&)      ;
    virtual QuantumNumber* Clone(        ) const;

    int Getm() const {return m;}

private:
    int m;
};

#endif //_Slater_m_QuantumNumer_H_
