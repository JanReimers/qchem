// File: Hamiltonian.C  Interface a Hamiltonianian operator.
export module qchem.Hamiltonian;
export import qchem.ChargeDensity;
import qchem.Streamable;
export import qchem.Energy;

//
//  Abstract base for any HamiltonianTerm (HT) terms in the Hamiltonian.
//  We have two distinct types of HT:
//  Static_HT - Does not depend on the Charge Dnesity (CD), and therefore does change during iterations
//  Dynamic_HT - Vee, Vac which depend on the CD, and change with each iteration.
//
export class Static_HT
    : public virtual Streamable
    , public virtual Static_CC
{
public:
    typedef SMatrix<double> SMat;
    typedef TOrbital_IBS<double> ibs_t;

    virtual const SMat& GetMatrix(const ibs_t*,const Spin&) const=0;
    virtual void        GetEnergy(EnergyBreakdown&,  const DM_CD*) const=0;
    virtual bool        IsPolarized() const {return false;}
};

export class Dynamic_HT
    : public virtual Streamable
    , public virtual Dynamic_CC
{
public:
    typedef SMatrix<double> SMat;
    typedef TOrbital_IBS<double> ibs_t;    
    virtual const SMat& GetMatrix(const ibs_t*,const Spin&,const DM_CD*) const=0; 
    virtual void        GetEnergy(EnergyBreakdown&,  const DM_CD*) const=0;
    virtual bool        IsPolarized() const {return false;}
};





export class Hamiltonian
    : public virtual Streamable
{
public:
    typedef SMatrix<double> SMat;
    typedef TOrbital_IBS<double> ibs_t;

    virtual void            Add             (      Static_HT*)      =0;
    virtual void            Add             (      Dynamic_HT*)      =0;
    virtual SMatrix<double>            GetMatrix(const ibs_t*,const Spin&,const DM_CD*)=0;
    virtual EnergyBreakdown GetTotalEnergy  (  const DM_CD*    ) const=0;
    virtual bool            IsPolarized() const=0;
};

