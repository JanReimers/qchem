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


struct DIISParams
{
    size_t Nproj;  //Number of terms to keep for proections.
    double EMax;   //DIIS starts when E<EMax
    double EMin;   //DIIS stops when E<EMin
    double SVTol;  //DIIS bails out when the minimum singular value of B matrix is < SVTol;
};

#include "oml/smatrix.h"
template <class T> class LinearSolver;
namespace oml
{
    template <class T> class SVDSolver;
}
class SCFIrrepAccelerator_DIIS : public virtual SCFIrrepAccelerator
{
public:
    SCFIrrepAccelerator_DIIS(const DIISParams&);
    virtual ~SCFIrrepAccelerator_DIIS();
    virtual void Init(const LASolver<double>*);
    virtual SMat Project(const SMat& F, const SMat& DPrime); //Fock and Density matrices.

private:
    typedef  Vector<double> RVec;

    SMat BuildB() const;
    bool IsSingular(const SMat& B, double SVtol) const;
    RVec SolveC(const SMat& B) const;
    SMat Project(const SMat& FPrime);
    void AppendAndPurge(const SMat& FPrime, const Mat& E, double En, size_t N);

    DIISParams                   itsParams; 
    std::vector<SMatrix<double>> itsFPrimes;
    std::vector< Matrix<double>> itsEs; //Error matrices
    std::vector<        double > itsEns; //Errors ||E||
    const LASolver   <double>*   itsLaSolver; //Knows the ortho transform
         LinearSolver<double>*   itsLinearSolver;
    oml::SVDSolver   <double>*   itsSVDSolver;

};


class SCFAccelerator_DIIS : public virtual SCFAccelerator
{
public:
    SCFAccelerator_DIIS(const DIISParams&);
    ~SCFAccelerator_DIIS();
    virtual SCFIrrepAccelerator* Create(const TOrbital_IBS<double>*) const;

private:
    DIISParams itsParams;
};

#endif //_SCFAccelerator_H_