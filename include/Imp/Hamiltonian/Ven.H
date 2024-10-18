// File: ExactVen.H  Nuclear potential.
#ifndef _ExactVen_H_
#define _ExactVen_H_




#include "Imp/Hamiltonian/HamiltonianTerm.H"
#include <Cluster.H>
#include "Misc/rc_ptr.H"

class ExactVen
    : public virtual HamiltonianTerm
    , private        HamiltonianTermImplementation
{
public:
    ExactVen(                         );
    ExactVen(const rc_ptr<Cluster>& cl);
    // Required by HamiltonianTerm
    virtual void GetEnergy       (TotalEnergy&               ) const;

    // Required by PMStreamableObject
    virtual std::ostream&   Write(std::ostream&) const;
    virtual std::istream&   Read (std::istream&)      ;

private:
    virtual bool DependsOnChargeDensity() const {return false;}
    virtual SMat CalculateHamiltonianMatrix(const IrrepBasisSet*,const Spin&) const;

    rc_ptr<Cluster> theCluster;
};

#endif //_ExactVen_H_
