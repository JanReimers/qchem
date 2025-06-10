// FIle: SCFAccelerator.H  Interface for an SCF accelerator alogrithm
#ifndef _SCFAccelerator_H_
#define _SCFAccelerator_H_

#include <types.H>
class SCFIrrepAccelerator
{
protected:
    typedef Vector<double>  RVec;
    typedef  Matrix<double>  Mat;
    typedef SMatrix<double> SMat;
public:
    virtual ~SCFIrrepAccelerator() {};
    virtual void Init(const LASolver<double>*)=0; //Solver knows how to transform to orthogonal basis.
    virtual void UseFD(const SMat& F, const SMat& DPrime)=0;
    virtual void SetProjection(const RVec& cs)=0;
    virtual  Mat CalculateError() =0;
    virtual SMat Project()=0; 
};

#include <iosfwd>
class SCFAccelerator
{
public:
    virtual ~SCFAccelerator() {};
    virtual SCFIrrepAccelerator* Create(const Irrep_QNs&)=0;
    virtual bool CalculateProjections()=0;
    virtual void ShowLabels     (std::ostream&) const=0;
    virtual void ShowConvergence(std::ostream&) const=0;
    virtual double GetError() const=0;
};


#endif //_SCFAccelerator_H_