// File: Atom/l/Slater_IEC.H Integral Engine Client for Slater atom basis functions.
#ifndef _Atom_l_Slater_IEC_H_
#define _Atom_l_Slater_IEC_H_

#include "Imp/BasisSet/AtomIEClient.H"

namespace Atoml
{
namespace Slater
{
    
struct IrrepIEClient
    : public virtual ::IrrepIEClient
    , public AtomIrrepIEClient
{
    IrrepIEClient() {};
    IrrepIEClient(size_t N) : AtomIrrepIEClient(N) {};
private:
    virtual double Norm(double e, size_t l) const;
};

}} //namespace Atom_l::Slater

#endif // _Atom_l_Slater_IEC_H_
