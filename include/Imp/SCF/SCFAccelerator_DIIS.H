// FIle: SCFAccelerator_DIIS.H  Direct Inversion of the Iterative Subspace (DIIS) algorithm
#ifndef _SCFAccelerator_DIIS_H_
#define _SCFAccelerator_DIIS_H_

#include "SCFAccelerator.H"

struct DIISParams
{
    size_t Nproj;  //Number of terms to keep for proections.
    double EMax;   //DIIS starts when E<EMax
    double EMin;   //DIIS stops when E<EMin
    double SVTol;  //DIIS bails out when the minimum singular value of B matrix is < SVTol;
};

#include <vector>

#include <Orbital_QNs.H>
#include "oml/smatrix.h"
#include "oml/vector.h"
class SCFAccelerator_DIIS;
class SCFIrrepAccelerator_DIIS : public virtual SCFIrrepAccelerator
{
public:
    typedef std::vector<Mat> mv_t; //matrix-vector type.
    typedef std::pair<SMat,double> md_t;
    SCFIrrepAccelerator_DIIS(const DIISParams&,const Irrep_QNs&);
    virtual ~SCFIrrepAccelerator_DIIS();
    
    virtual void Init(const LASolver<double>*);
    virtual void UseFD(const SMat& F, const SMat& DPrime);
    virtual void SetProjection(const RVec& cs) {itsCs=cs;}
    virtual  Mat CalculateError() ;
    virtual SMat Project(); 
private:
    friend class SCFAccelerator_DIIS;
    


    double GetError() const {return itsLastEn;}
    // double GetSVMin() const {return itsLastSVMin;}
    size_t GetNproj() const {return itsEs.size();}
    // bool   Bailout () const {return itsBailout;}
    void   GlobalBailout() {itsBailout=true;}
    double GetError(size_t i, size_t j) const {return Dot(itsEs[i],itsEs[j]);}

    void Append(const SMat& FPrime, const Mat& E, double En);
    void AppendFPrime() {Append(itsFPrime,itsLastE,itsLastEn);}
    void Purge1();
    // RVec Solve();


    typedef  Vector<double> RVec;

    
    DIISParams                   itsParams; 
    Irrep_QNs                    itsIrrep;
    double                       itsLastEn;  // [F',D']
    // double                       itsLastSVMin;  //Smallest singular value for latest B matric.
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
    virtual SCFIrrepAccelerator* Create(const Irrep_QNs&);
    virtual bool CalculateProjections();
    virtual void ShowLabels     (std::ostream&) const;
    virtual void ShowConvergence(std::ostream&) const;
    virtual double GetError() const;

private:
    static double GetMinSV(const SMatrix<double>& B);
    static RVec   SolveC(const SMatrix<double>& B);
    SMat BuildB() const;
    md_t BuildPrunedB(double svmin);
    void Purge1();
    bool BailoutChildren();
    size_t GetNProj() const;

  
    DIISParams itsParams;
    std::vector<SCFIrrepAccelerator_DIIS*> itsIrreps;

    bool itsBailout;
    double itsEn,itsLastSVMin;
    RVec itsCs;
};

#endif //_SCFAccelerator_DIIS_H_