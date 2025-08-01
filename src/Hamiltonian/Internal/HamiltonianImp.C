// File: HamiltonianImp.C  General matrix implementation of a Hamiltonian operator.
module;
#include <vector>
#include <memory>
export module qchem.Hamiltonian.Internal.Hamiltonian;
export import qchem.Hamiltonian;
import qchem.Cluster;

export class HamiltonianImp
    : public virtual Hamiltonian
{
public:
    HamiltonianImp();
    virtual void Add(Static_HT* );
    virtual void Add(Dynamic_HT*);

    virtual SMatrix<double>            GetMatrix(const ibs_t*,const Spin& S,const DM_CD*);
    virtual EnergyBreakdown GetTotalEnergy  (const DM_CD* ) const;
    virtual bool            IsPolarized() const {return itsIsPolarized;}
    virtual std::ostream&   Write(std::ostream&) const;
 
protected:
    typedef std::shared_ptr<const Cluster> cl_t;
    void InsertStandardTerms(const cl_t & cl);
    typedef std::vector<std::unique_ptr< Static_HT>> shtv_t;
    typedef std::vector<std::unique_ptr<Dynamic_HT>> dhtv_t;

    shtv_t itsSHTs;
    dhtv_t itsDHTs;
    bool   itsIsPolarized;
};
