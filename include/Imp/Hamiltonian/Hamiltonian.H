// File: HamiltonianImplementation.H  General matrix implementation of a Hamiltonian operator.
#ifndef _HamiltonianImp_H_
#define _HamiltonianImp_H_



#include <Hamiltonian.H>
#include "Imp/Containers/ptr_vector.h"
#include <memory>

class HamiltonianImp
    : public virtual Hamiltonian
{
public:
    HamiltonianImp();
    virtual void Add(Static_HT* );
    virtual void Add(Dynamic_HT*);

    virtual SMat        GetMatrix(const ibs_t*,const Spin& S,const DM_CD*);
    virtual EnergyBreakdown GetTotalEnergy  (const DM_CD* ) const;
    
    virtual std::ostream&     Write(std::ostream&) const;
 
protected:
    typedef std::shared_ptr<const Cluster> cl_t;
    void InsertStandardTerms(const cl_t & cl);
    optr_vector1< Static_HT*> itsSHTs;
    optr_vector1<Dynamic_HT*> itsDHTs;
};

#endif //_HamiltonianImp_H_
