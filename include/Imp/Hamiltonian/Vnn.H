// File: Vnn.H  Nuclear-Nuclear potential.
#ifndef _Vnn_H_
#define _Vnn_H_



#include "Imp/Hamiltonian/HamiltonianTerm.H"
#include <Cluster.H>
#include <memory>

class Vnn
    : public virtual HamiltonianTerm
    , private        HamiltonianTermImp
{
public:
    typedef std::shared_ptr<const Cluster> cl_t;
    Vnn(                         );
    Vnn(cl_t& cl);
    // Required by HamiltonianTerm
    virtual void GetEnergy       (TotalEnergy&               ) const;
    // Required by PMStreamableObject
    virtual std::ostream& Write(std::ostream&) const;
    virtual std::istream& Read (std::istream&)      ;

private:
    virtual bool DependsOnChargeDensity() const {return false;}
    virtual SMat CalculateHamiltonianMatrix(const IrrepBasisSet*,const Spin&) const;
    cl_t theCluster;
};

#endif //_Vnn_H_
