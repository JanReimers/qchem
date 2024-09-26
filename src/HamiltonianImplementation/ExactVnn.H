// File: ExactVnn.H  Nuclear potential.
#ifndef _ExactVnn_H_
#define _ExactVnn_H_



#include "Cluster/Cluster.H"
#include "Misc/rc_ptr.H"
#include "HamiltonianImplementation/HamiltonianTermImplementation.H"

class ExactVnn
    : public virtual HamiltonianTerm
    , private        HamiltonianTermImplementation
{
public:
    ExactVnn(                         );
    ExactVnn(const rc_ptr<Cluster>& cl);
    // Required by HamiltonianTerm
    virtual void GetEnergy       (TotalEnergy&               ) const;
    // Required by PMStreamableObject
    virtual std::ostream& Write(std::ostream&) const;
    virtual std::istream& Read (std::istream&)      ;

private:
    virtual bool DependsOnChargeDensity() const {return false;}
    virtual SMat CalculateHamiltonianMatrix(const BasisSet*,const Spin&) const;
    rc_ptr<Cluster> theCluster;
};

#endif //_ExactVnn_H_
