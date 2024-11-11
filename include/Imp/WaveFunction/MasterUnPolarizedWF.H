// File: MasterUnPolarizedWF.H  Wave function for an unpolarized atom.
#ifndef _MasterUnPolarizedWF_H_
#define _MasterUnPolarizedWF_H_



#include "WaveFunction.H"

class MasterUnPolarizedWF
    : public virtual WaveFunction
{
public:
    MasterUnPolarizedWF();
    MasterUnPolarizedWF(const BasisSet*,const ElectronConfiguration*);
    virtual ~MasterUnPolarizedWF();

    virtual void           DoSCFIteration      (Hamiltonian&   )      ;
    virtual ChargeDensity* GetChargeDensity    (Spin           ) const;
    virtual void           UpdateElectronDumper(ElectronDumper&)      ;
    virtual void           FillOrbitals        (const ElectronConfiguration*, const Spin&);
    virtual SCFIterator*   MakeIterator        (Hamiltonian*, ChargeDensity* guess, double nElectrons);
    virtual void           DisplayEigen        (               ) const;


    virtual std::ostream&       Write(std::ostream&) const;
    virtual std::istream&       Read (std::istream&)      ;

private:
    WaveFunction* itsGroup;
    const ElectronConfiguration* itsEC;

};

#endif //_MasterUnPolarizedWF_H_
