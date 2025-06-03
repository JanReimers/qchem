// File: SCFIterator.H  Interface for an object that manages SCF convergence.


#include "types.H"

class SCFIterator
{
public:
    SCFIterator(const BasisSet*, const ElectronConfiguration*, Hamiltonian*,DM_CD* cd=0);
    virtual ~SCFIterator();
    virtual bool Iterate(const SCFIterationParams& ipar);

    const WaveFunction* GetWaveFunction() const {return itsWaveFunction;}
    EnergyBreakdown     GetEnergy() const;
    DM_CD*              GetChargeDensity () const;

private:
    void Initialize(DM_CD* cd);  //Does on iteration to set up the exact charge density.
    void DisplayEnergies(int i, const EnergyBreakdown&,  double relax, double chargeDensityChange, double fitError) const;
    void DisplayEigen   () const;

    //All owned, see destructor.
    Hamiltonian*  itsHamiltonian;
    WaveFunction* itsWaveFunction;  
    DM_CD*        itsCD;
    DM_CD*        itsOldCD;
};


#include <vector>
class SCFIrrepIterator_DIIS
{
public:
    SCFIrrepIterator_DIIS(Orbitals* os,const SMatrix<double>&);

private:
    std::vector<SMatrix<double>> itsF_Primes;


};

#include <memory>
//  Accelerated convergence algo described in: https://www.chem.fsu.edu/~deprince/programming_projects/diis/
class SCFIterator_DIIS : public SCFIterator
{
public:
    SCFIterator_DIIS(const BasisSet* bs, const ElectronConfiguration* ec, Hamiltonian*,DM_CD* cd=0);
    ~SCFIterator_DIIS();
    // virtual bool Iterate(const SCFIterationParams& ipar);

private:
    typedef std::shared_ptr<SCFIrrepIterator_DIIS> ii_t; //Irrep iterator type
    std::vector<ii_t> itsSCFIterators;
    
};

