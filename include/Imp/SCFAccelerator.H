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
    virtual void Init(const LASolver<double>*)=0; //Solver knows how to transform to orthogonal basis.
    virtual SMat Project(const SMat& F, const SMat& DPrime)=0; //Fock and Density matrices.

};

class SCFAccelerator
{
public:
    virtual ~SCFAccelerator() {};
    virtual SCFIrrepAccelerator* Create(const TOrbital_IBS<double>*) const=0;
};


#include <deque>
#include "oml/smatrix.h"

class SCFIrrepAccelerator_DIIS : public virtual SCFIrrepAccelerator
{
public:
    SCFIrrepAccelerator_DIIS(const SMat& overlap);
    virtual ~SCFIrrepAccelerator_DIIS();
    virtual void Init(const LASolver<double>*);
    virtual SMat Project(const SMat& F, const SMat& DPrime); //Fock and Density matrices.

private:
    typedef  Vector<double> RVec;

    SMat BuildB() const;
    bool IsSingular(const SMat& B, double SVtol) const;
    RVec SolveC(const SMat& B) const;
    SMat Solve(const SMat& Fprime);
    void AppendAndPurge(const SMat& FPrime, const Mat& E, double En, size_t N);
    
    std::deque<SMatrix<double>> itsFPrimes;
    std::deque< Matrix<double>> itsEs; //Error matrices
    std::vector<       double > itsEns; //Errors ||E||
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