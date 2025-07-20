// FIle: SCFAccelerator.C  Interface for an SCF accelerator alogrithm
module;
#include <iosfwd>
export module qchem.SCFAccelerator;
export import oml;
export import qchem.Symmetry.Irrep;
export import qchem.LASolver;

export class SCFIrrepAccelerator
{
protected:
    typedef SMatrix<double> SMat;
public:
    virtual ~SCFIrrepAccelerator() {};
    virtual void UseFD(const SMat& F, const SMat& DPrime)=0;
    virtual SMat Project()=0; 
};

export class SCFAccelerator
{
public:
    virtual ~SCFAccelerator() {};
    virtual SCFIrrepAccelerator* Create(const LASolver<double>*,const Irrep_QNs&, int occ)=0;
    virtual bool CalculateProjections()=0;
    virtual void ShowLabels     (std::ostream&) const=0;
    virtual void ShowConvergence(std::ostream&) const=0;
    virtual double GetError() const=0;
};


