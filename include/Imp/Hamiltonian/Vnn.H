// File: Vnn.H  Nuclear-Nuclear potential.
#ifndef _Vnn_H_
#define _Vnn_H_



#include "Imp/Hamiltonian/HamiltonianTerm.H"
#include <Cluster.H>
#include "Imp/Misc/rc_ptr.H"

class Vnn
    : public virtual HamiltonianTerm
    , private        HamiltonianTermImp
{
public:
    Vnn(                         );
    Vnn(const rc_ptr<Cluster>& cl);
    // Required by HamiltonianTerm
    virtual void GetEnergy       (TotalEnergy&               ) const;
    // Required by PMStreamableObject
    virtual std::ostream& Write(std::ostream&) const;
    virtual std::istream& Read (std::istream&)      ;

private:
    virtual bool DependsOnChargeDensity() const {return false;}
    virtual SMat CalculateHamiltonianMatrix(const IrrepBasisSet*,const Spin&) const;
    rc_ptr<Cluster> theCluster;
};

#endif //_Vnn_H_
