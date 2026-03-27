// FIle: SCFAcceleratorDIIS.H  Direct Inversion of the Iterative Subspace (DIIS) algorithm
module;
#include <deque>
#include <vector>
#include <iostream>
#include <blaze/Math.h>
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
    typedef std::deque< rmat_t> mv_t; //matrix-vector type.
    typedef std::deque<rsmat_t> sv_t; //smatrix-vector type.
    typedef std::deque<double> dv_t ; //doubles
    
    SCFIrrepAcceleratorDIIS(const DIISParams&,const LASolver_blaze<double>*,const Irrep_QNs&,const rvec_t& cs);
    virtual ~SCFIrrepAcceleratorDIIS();
    
    virtual void UseFD(const rsmat_t& F, const rsmat_t& DPrime);
    virtual rsmat_t Project(); 
private:
    friend class SCFAcceleratorDIIS;
    size_t GetNproj() const {return itsEs.size();}
    double GetError() const {return itsEn;}
    double GetError(size_t i, size_t j) const {return sum(itsEs[i] % itsEs[j]);}
    void Append1();
    void Purge1();
    
    DIISParams itsParams; 
    Irrep_QNs  itsIrrep;
    // All of these are the the most recent values
    rsmat_t    itsFPrime,itsDPrime;    
    rmat_t     itsE;
    double     itsEn;  // [F',D']
    // Caches for F',E,En
    sv_t itsFPrimes;
    mv_t itsEs; //Error matrices [F',D']
    dv_t itsEns; //Errors ||E||=FNorm[F',D']
    

    const rvec_t& itsCs;  //Projection coefficients from SCFAcceleratorDIIS class.

    const LASolver_blaze   <double>*   itsLaSolver_blaze; //Knows the ortho transform

};


export class SCFAcceleratorDIIS : public virtual SCFAccelerator
{
public:

    SCFAcceleratorDIIS(const DIISParams&);
    ~SCFAcceleratorDIIS();
    virtual SCFIrrepAccelerator* Create(const LASolver_blaze<double>*,const Irrep_QNs&, int occ);
    virtual bool   CalculateProjections();
    virtual void   ShowLabels     (std::ostream&) const;
    virtual void   ShowConvergence(std::ostream&) const;
    virtual double GetError() const;

private:
    typedef std::vector<rmat_t> mv_t; //matrix-vector type.


    static double GetMinSV(const rsmat_t& B);
    static rvec_t SolveC  (const rsmat_t& B);
    
    struct md_t{rsmat_t B;double sv;};

    md_t    BuildB() const;
    rsmat_t BuildPrunedB(double svmin);
    size_t  Append1();
    size_t  Purge1();
    size_t  GetNProj() const;
    bool    HasProjection() const {return itsCs.size()>=2;}

  
    DIISParams itsParams;
    std::vector<SCFIrrepAcceleratorDIIS*> itsIrreps;

    double itsEn,itsLastSVMin;
    rvec_t itsCs;
};

