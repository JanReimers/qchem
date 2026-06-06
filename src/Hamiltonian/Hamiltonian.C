// File: Hamiltonian.C  Interface a Hamiltonianian operator.
export module qchem.Hamiltonian;
export import qchem.ChargeDensity;
import qchem.Streamable;
export import qchem.Energy;
export import qchem.Hamiltonian.Types;


export namespace qchem::Hamiltonian
{

using ChargeDensity::DM_CD;

//
//  Abstract base for any HamiltonianTerm (HT) terms in the Hamiltonian.
//  We have two distinct types of HT:
//  Static_HT - Does not depend on the Charge Dnesity (CD), and therefore does change during iterations
//  Dynamic_HT - Vee, Vac which depend on the CD, and change with each iteration.
//
class Static_HT
    : public virtual Streamable
    , public virtual ChargeDensity::Static_CC
{
public:
    virtual const rsmat_t& GetMatrix(const obs_t*,const Spin&) const=0;
    virtual void           GetEnergy(EnergyBreakdown&,  const DM_CD*) const=0;
    virtual bool           IsPolarized   () const {return false;}
    virtual bool           IsRelativistic() const {return false;}
};

class Dynamic_HT
    : public virtual Streamable
    , public virtual ChargeDensity::Dynamic_CC
{
public:
    virtual const rsmat_t& GetMatrix(const obs_t*,const Spin&,const DM_CD*) const=0; 
    virtual void           GetEnergy(EnergyBreakdown&,  const DM_CD*) const=0;
    virtual bool           IsPolarized   () const {return false;}
    virtual bool           IsRelativistic() const {return false;}
};





class Hamiltonian
    : public virtual Streamable
{
public:
    virtual void            Add             ( Static_HT*)=0;
    virtual void            Add             (Dynamic_HT*)=0;
    virtual rsmat_t         GetMatrix(const obs_t*,const Spin&,const DM_CD*)=0;
    virtual EnergyBreakdown GetTotalEnergy  (  const DM_CD*    ) const=0;
    virtual bool            IsPolarized   () const=0;
    virtual bool            IsRelativistic() const=0;
};

} //namespace

