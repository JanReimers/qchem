// FIle: SCFAcceleratorDIIS.H  Direct Inversion of the Iterative Subspace (DIIS) algorithm
module;
#include <deque>
#include <vector>
#include <iostream>
export module qchem.SCFAccelerator.Internal.SCFAcceleratorDIIS;
export import qchem.SCFAccelerator;

export struct DIISParams
{
    size_t Nproj;  //Number of terms to keep for proections.
    double EMax;   //DIIS starts when E<EMax
    double EMin;   //DIIS stops when E<EMin
    double SVTol;  //DIIS bails out when the minimum singular value of B matrix is < SVTol;
};



export class SCFAcceleratorDIIS;
class SCFIrrepAcceleratorDIIS : public virtual SCFIrrepAccelerator
{
public:
    typedef Matrix<double> Mat;
    typedef std::deque<Mat   > mv_t; //matrix-vector type.
    typedef std::deque< SMatrix<double>   > sv_t; //smatrix-vector type.
    typedef std::deque<double> dv_t ; //doubles
    
    SCFIrrepAcceleratorDIIS(const DIISParams&,const LASolver<double>*,const Irrep_QNs&,const RVec& cs);
    virtual ~SCFIrrepAcceleratorDIIS();
    
    virtual void UseFD(const SMatrix<double>& F, const SMatrix<double>& DPrime);
    virtual SMatrix<double> Project(); 
private:
    friend class SCFAcceleratorDIIS;
    size_t GetNproj() const {return itsEs.size();}
    double GetError() const {return itsEn;}
    double GetError(size_t i, size_t j) const {return Dot(itsEs[i],itsEs[j]);}
    void Append1();
    void Purge1();
    
    DIISParams itsParams; 
    Irrep_QNs  itsIrrep;
    // All of these are the the most recent values
    SMatrix<double>   itsFPrime,itsDPrime;    
    Mat    itsE;
    double itsEn;  // [F',D']
    // Caches for F',E,En
    sv_t itsFPrimes;
    mv_t itsEs; //Error matrices [F',D']
    dv_t itsEns; //Errors ||E||=FNorm[F',D']
    

    const RVec&                  itsCs;  //Projection coefficients from SCFAcceleratorDIIS class.

    const LASolver   <double>*   itsLaSolver; //Knows the ortho transform

};


export class SCFAcceleratorDIIS : public virtual SCFAccelerator
{
public:

    SCFAcceleratorDIIS(const DIISParams&);
    ~SCFAcceleratorDIIS();
    virtual SCFIrrepAccelerator* Create(const LASolver<double>*,const Irrep_QNs&, int occ);
    virtual bool   CalculateProjections();
    virtual void   ShowLabels     (std::ostream&) const;
    virtual void   ShowConvergence(std::ostream&) const;
    virtual double GetError() const;

private:
    typedef  Matrix<double>  Mat;
    typedef SMatrix<double> SMat;
    typedef std::vector< Mat> mv_t; //matrix-vector type.


    static double GetMinSV(const SMat& B);
    static RVec   SolveC  (const SMat& B);
    
    struct md_t{SMat B;double sv;};

    md_t   BuildB() const;
    SMatrix<double>   BuildPrunedB(double svmin);
    size_t Append1();
    size_t Purge1();
    size_t GetNProj() const;
    bool   HasProjection() const {return itsCs.size()>=2;}

  
    DIISParams itsParams;
    std::vector<SCFIrrepAcceleratorDIIS*> itsIrreps;

    double itsEn,itsLastSVMin;
    RVec   itsCs;
};

