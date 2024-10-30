// File: Slater_m/QuantumNumer.H  Spherical harmonic Ylm symmetry.
#ifndef _Slater_m_QuantumNumer_H_
#define _Slater_m_QuantumNumer_H_



#include <QuantumNumber.H>

//---------------------------------------------------------------------------------
//
//  angular momentum L and m.
//

class YlmQN
    : public QuantumNumber
{
public:
    YlmQN(            );
    YlmQN(int l, int m);

    virtual bool Match(const QuantumNumber&) const;
    virtual int GetDegeneracy() const;

    virtual std::ostream&       Write(std::ostream&) const;
    virtual std::istream&       Read (std::istream&)      ;
    virtual QuantumNumber* Clone(        ) const;

    int Getl() const {return l;}
    int Getm() const {return m;}

private:
    int l,m;
};

#endif //_Slater_m_QuantumNumer_H_
