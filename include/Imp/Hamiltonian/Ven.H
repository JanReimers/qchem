// File: Ven.H  Electron-Nuclear potential.
#ifndef _Ven_H_
#define _Ven_H_




#include "Imp/Hamiltonian/HamiltonianTerm.H"
#include <Cluster.H>
#include <memory>

class Ven
    : public virtual Static_HT
    , private        Static_HT_Imp
{
public:
    typedef std::shared_ptr<const Cluster> cl_t;
    Ven(                         );
    Ven(const cl_t& cl);
    // Required by HamiltonianTerm
    virtual void GetEnergy       (TotalEnergy&,const DM_CD* cd               ) const;

    // Required by PMStreamableObject
    virtual std::ostream&   Write(std::ostream&) const;
    virtual std::istream&   Read (std::istream&)      ;

private:
    virtual SMat CalculateMatrix(const ibs_t*,const Spin&) const;

    cl_t theCluster;
};

#endif //_Ven_H_
