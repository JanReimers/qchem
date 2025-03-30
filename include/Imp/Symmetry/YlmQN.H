// File: Symmetry/YlmQN.H  Spherical harmonic Ylm symmetry.
#ifndef _Symmetry_YlmQN_H_
#define _Symmetry_YlmQN_H_

#include "Imp/Symmetry/YlQN.H"


//---------------------------------------------------------------------------------
//
//  angular momentum L and m.
//

class YlmQN
    : public virtual QNs
    , public YlQN
{
public:
    YlmQN(            );
    YlmQN(int l, int m);

    virtual bool   Match(const QNs&) const;
    virtual int    GetDegeneracy() const;
    QNs* AddPrincipleQN(int index) const;
    virtual std::pair<int,int> GetN(const int (&N)[4], const int (&Nv)[4], int NUnpaired) const;

    virtual std::ostream&  Write(std::ostream&) const;
    virtual std::istream&  Read (std::istream&)      ;
    virtual AngularQN* Clone(        ) const;

    int Getm() const {return m;}
    using YlQN::GetL;
protected:
    int m;
};

#endif //_Symmetry_YlmQN_H_
