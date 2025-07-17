// FIle: SCFAccelerator_DIIS.H  Direct Inversion of the Iterative Subspace (DIIS) algorithm
#ifndef _SCFAccelerator_DIIS_H_
#define _SCFAccelerator_DIIS_H_

#include <deque>
#include <vector>
#include <SCFAccelerator/SCFAccelerator.H>
#include <Symmetry/Irrep_QNs.H>
import oml;

struct DIISParams
{
    size_t Nproj;  //Number of terms to keep for proections.
    double EMax;   //DIIS starts when E<EMax
    double EMin;   //DIIS stops when E<EMin
    double SVTol;  //DIIS bails out when the minimum singular value of B matrix is < SVTol;
};



class SCFAccelerator_DIIS;
class SCFIrrepAccelerator_DIIS : public virtual SCFIrrepAccelerator
{
public:
    typedef Matrix<double> Mat;
    typedef Vector<double> RVec;
    typedef std::deque<Mat   > mv_t; //matrix-vector type.
    typedef std::deque<SMat  > sv_t; //smatrix-vector type.
    typedef std::deque<double> dv_t ; //doubles
    
    SCFIrrepAccelerator_DIIS(const DIISParams&,const LASolver<double>*,const Irrep_QNs&,const RVec& cs);
    virtual ~SCFIrrepAccelerator_DIIS();
    
    virtual void UseFD(const SMat& F, const SMat& DPrime);
    virtual SMat Project(); 
private:
    friend class SCFAccelerator_DIIS;
    size_t GetNproj() const {return itsEs.size();}
    double GetError() const {return itsEn;}
    double GetError(size_t i, size_t j) const {return Dot(itsEs[i],itsEs[j]);}
    void Append1();
    void Purge1();
    
    DIISParams itsParams; 
    Irrep_QNs  itsIrrep;
    // All of these are the the most recent values
    SMat   itsFPrime,itsDPrime;    
    Mat    itsE;
    double itsEn;  // [F',D']
    // Caches for F',E,En
    sv_t itsFPrimes;
    mv_t itsEs; //Error matrices [F',D']
    dv_t itsEns; //Errors ||E||=FNorm[F',D']
    

    const RVec&                  itsCs;  //Projection coefficients from SCFAccelerator_DIIS class.

    const LASolver   <double>*   itsLaSolver; //Knows the ortho transform

};


class SCFAccelerator_DIIS : public virtual SCFAccelerator
{
public:

    SCFAccelerator_DIIS(const DIISParams&);
    ~SCFAccelerator_DIIS();
    virtual SCFIrrepAccelerator* Create(const LASolver<double>*,const Irrep_QNs&, int occ);
    virtual bool   CalculateProjections();
    virtual void   ShowLabels     (std::ostream&) const;
    virtual void   ShowConvergence(std::ostream&) const;
    virtual double GetError() const;

private:
    typedef Vector<double>  RVec;
    typedef  Matrix<double>  Mat;
    typedef SMatrix<double> SMat;
    typedef std::vector< Mat> mv_t; //matrix-vector type.


    static double GetMinSV(const SMat& B);
    static RVec   SolveC  (const SMat& B);
    
    struct md_t{SMat B;double sv;};

    md_t   BuildB() const;
    SMat   BuildPrunedB(double svmin);
    size_t Append1();
    size_t Purge1();
    size_t GetNProj() const;
    bool   HasProjection() const {return itsCs.size()>=2;}

  
    DIISParams itsParams;
    std::vector<SCFIrrepAccelerator_DIIS*> itsIrreps;

    double itsEn,itsLastSVMin;
    RVec   itsCs;
};

#endif //_SCFAccelerator_DIIS_H_