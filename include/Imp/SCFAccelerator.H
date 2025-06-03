// FIle: SCFAccelerator.H  Interface for an accelerator alogrithm
#ifndef _SCFAccelerator_H_
#define _SCFAccelerator_H_

#include <types.H>
class SCFIrrepAccelerator
{
protected:
    typedef  Matrix<double>  Mat;
    typedef SMatrix<double> SMat;
public:
    virtual ~SCFIrrepAccelerator() {};
    virtual void Init(const LASolver<double>*)=0;;

};

class SCFAccelerator
{
public:
    virtual ~SCFAccelerator() {};
    virtual SCFIrrepAccelerator* Create(const TOrbital_IBS<double>*) const=0;
};


#include <vector>
#include "oml/smatrix.h"

class SCFIrrepAccelerator_DIIS : public virtual SCFIrrepAccelerator
{
public:
    SCFIrrepAccelerator_DIIS(const SMat& overlap);
    virtual ~SCFIrrepAccelerator_DIIS();
    virtual void Init(const LASolver<double>*);
private:
    std::vector<SMatrix<double>> itsF_Primes;
    SMat S;
    const LASolver<double>* itsLaSolver; //Knows the ortho transform

};


class SCFAccelerator_DIIS : public virtual SCFAccelerator
{
public:
    SCFAccelerator_DIIS();
    ~SCFAccelerator_DIIS();
    virtual SCFIrrepAccelerator* Create(const TOrbital_IBS<double>*) const;

private:

};

#endif //_SCFAccelerator_H_