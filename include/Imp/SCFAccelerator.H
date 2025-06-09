// FIle: SCFAccelerator.H  Interface for an accelerator alogrithm
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
    virtual void Init(const LASolver<double>*, const Irrep_QNs&)=0; //Solver knows how to transform to orthogonal basis.
    virtual void UseFD(const SMat& F, const SMat& DPrime)=0;
    virtual void SetProjection(const RVec& cs)=0;
    virtual  Mat CalculateError() =0;
    virtual SMat Project(const SMat& F, const SMat& DPrime)=0; //Fock and Density matrices.
    virtual SMat Project()=0; 
};

#include <iosfwd>
class SCFAccelerator
{
public:
    virtual ~SCFAccelerator() {};
    virtual SCFIrrepAccelerator* Create(const TOrbital_IBS<double>*)=0;
    virtual bool CalculateProjections()=0;
    virtual void ShowLabels     (std::ostream&) const=0;
    virtual void ShowConvergence(std::ostream&) const=0;
    virtual double GetError() const=0;
};


struct DIISParams
{
    enum {global,irrep} type;
    size_t Nproj;  //Number of terms to keep for proections.
    double EMax;   //DIIS starts when E<EMax
    double EMin;   //DIIS stops when E<EMin
    double SVTol;  //DIIS bails out when the minimum singular value of B matrix is < SVTol;
};

#include <vector>

#include <Orbital_QNs.H>
#include "oml/smatrix.h"
#include "oml/vector.h"
template <class T> class LinearSolver;
namespace oml
{
    template <class T> class SVDSolver;
}
class SCFIrrepAccelerator_DIIS : public virtual SCFIrrepAccelerator
{
public:
    typedef std::vector<Mat> mv_t; //matrix-vector type.
    typedef std::pair<SMat,double> md_t;
    SCFIrrepAccelerator_DIIS(const DIISParams&);
    virtual ~SCFIrrepAccelerator_DIIS();
    virtual void Init(const LASolver<double>*, const Irrep_QNs&);
    virtual void UseFD(const SMat& F, const SMat& DPrime);
    virtual void SetProjection(const RVec& cs) {itsCs=cs;}

    virtual  Mat CalculateError() ;
    virtual bool CalculateProjections();
    virtual SMat Project(const SMat& F, const SMat& DPrime); //Fock and Density matrices.
    virtual SMat Project(); 
   
    double GetError() const {return itsLastEn;}
    double GetSVMin() const {return itsLastSVMin;}
    size_t GetNproj() const {return itsEs.size();}
    bool   Bailout () const {return itsBailout;}
    void   GlobalBailout() {itsBailout=true;}
    void   AppendFPrime() {Append(itsFPrime,itsLastE,itsLastEn);}
    SMat   BuildRawB() const {return BuildRawB(itsEs);}


    static SMat&  AddBEdges(SMat& B);
    static SMat   BuildRawB(const mv_t&);
    static SMat   BuildB(const mv_t&);
    static double GetMinSV(const SMatrix<double>& B);
    static RVec   SolveC(const SMatrix<double>& B);
    void Purge1();

private:
    typedef  Vector<double> RVec;

    md_t BuildPrunedB(const mv_t&,double svmin);
    // RVec SolveC(const SMat& B) const;
    RVec Solve();
    SMat Project(const RVec& c);
    void Append(const SMat& FPrime, const Mat& E, double En);

    DIISParams                   itsParams; 
    Irrep_QNs                    itsIrrep;
    double                       itsLastEn;  // [F',D']
    double                       itsLastSVMin;  //Smallest singular value for latest B matric.
    SMat                         itsFPrime,itsDPrime;    
    Mat                          itsLastE;
    RVec                         itsCs;  //Projection coefficients.
    bool                         itsBailout;

    std::vector<SMatrix<double>> itsFPrimes;
    mv_t                         itsEs; //Error matrices
    std::vector<        double > itsEns; //Errors ||E||
    const LASolver   <double>*   itsLaSolver; //Knows the ortho transform

};


class SCFAccelerator_DIIS : public virtual SCFAccelerator
{
protected:
    typedef Vector<double>  RVec;
    typedef  Matrix<double>  Mat;
    typedef SMatrix<double> SMat;
public:
    typedef std::vector< Mat> mv_t; //matrix-vector type.
    typedef std::pair<SMat,double> md_t;

    SCFAccelerator_DIIS(const DIISParams&);
    ~SCFAccelerator_DIIS();
    virtual SCFIrrepAccelerator* Create(const TOrbital_IBS<double>*);
    virtual bool CalculateProjections();
    virtual void ShowLabels     (std::ostream&) const;
    virtual void ShowConvergence(std::ostream&) const;
    virtual double GetError() const;

private:
    md_t BuildPrunedB(double svmin);
    void Purge1();
    bool BailoutChildren();

  
    DIISParams itsParams;
    std::vector<SCFIrrepAccelerator_DIIS*> itsIrreps;

    bool itsBailout;
    double itsEn,itsLastSVMin;
    size_t itsN;
    RVec itsCs;
};

#endif //_SCFAccelerator_H_