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

    virtual void        Add             (      Static_HT* )      ;
    virtual void        Add             (      Dynamic_HT*)      ;
    virtual SMat        GetMatrix(const ibs_t*,const Spin& S,const DM_CD*);
    virtual TotalEnergy GetTotalEnergy  (const DM_CD* ) const;
    
    virtual std::ostream&     Write(std::ostream&) const;
    virtual std::istream&     Read (std::istream&)      ;

protected:
    void        UseChargeDensity(const DM_CD*  )      ;
    typedef std::shared_ptr<const Cluster> cl_t;
    void InsertStandardTerms(cl_t & cl);
    optr_vector1<Static_HT*> itsSHTs;
    optr_vector1<Dynamic_HT*> itsDHTs;
    const DM_CD*    itsCD;    //Exact charge density.
};

#endif //_HamiltonianImp_H_
