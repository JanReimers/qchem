// FIle: SCFAccelerator.H  Interface for an SCF accelerator alogrithm
#ifndef _SCFAccelerator_H_
#define _SCFAccelerator_H_

#include <types.H>
class SCFIrrepAccelerator
{
protected:
    typedef SMatrix<double> SMat;
public:
    virtual ~SCFIrrepAccelerator() {};
    virtual void UseFD(const SMat& F, const SMat& DPrime)=0;
    virtual SMat Project()=0; 
};

#include <iosfwd>
class SCFAccelerator
{
public:
    virtual ~SCFAccelerator() {};
    virtual SCFIrrepAccelerator* Create(const LASolver<double>*,const Irrep_QNs&, int occ)=0;
    virtual bool CalculateProjections()=0;
    virtual void ShowLabels     (std::ostream&) const=0;
    virtual void ShowConvergence(std::ostream&) const=0;
    virtual double GetError() const=0;
};


#endif //_SCFAccelerator_H_