// File: NumricalIEImp.H  Numerical integral engines.
#ifndef _NumricalIEImp_H_
#define _NumricalIEImp_H_



#include "NumericalIE.H"
#include "Misc/Persistent/IDRef.H"
#include "Mesh/MeshIntegrator.H"
#include "oml/vector.h"

template <class T> class NumericalIEImp
    : public virtual NumericalIE<T>
{
    typedef NumericalIE<T> Base;
    typedef typename Base::Vec   Vec;
    typedef typename Base::Mat   Mat;
    typedef typename Base::RVec  RVec;
    typedef typename Base::SMat  SMat;
    typedef typename Base::bs_t bs_t;
public:
    NumericalIEImp(                  );
    NumericalIEImp(const NumericalIEImp&);
    NumericalIEImp(Mesh*             );
    ~NumericalIEImp(                  );

    virtual void Insert(bs_t*);

    virtual SMat MakeOverlap    (                             ) const;
    virtual Mat  MakeOverlap    (bs_t& a,bs_t& b          ) const;
    virtual Vec  MakeOverlap    (const ScalarFunction<double>&) const;

    virtual SMat MakeRepulsion  (                             ) const;
    virtual Mat  MakeRepulsion  (bs_t& a,bs_t& b    ) const;
    virtual Vec  MakeRepulsion  (const ScalarFunction<double>&) const;

    virtual RVec MakeNormalization(                 ) const;
    virtual RVec MakeCharge       (                 ) const;
#ifdef USE_FOR_DEBUGGING_ANALYTIC
    virtual SMat MakeKinetic  (                     ) const;
    virtual SMat MakeNuclear  (const Cluster&       ) const;
    virtual void MakeOverlap3C  (ERI3&,bs_t&   ) const;
    virtual void MakeRepulsion3C(ERI3&,bs_t&   ) const;
    virtual void MakeRepulsion4C(ERI4&,ERIList&,const BasisGroup*   ) const;
#endif

    virtual std::ostream&           Write(std::ostream&) const;
    virtual std::istream&           Read (std::istream&)      ;
    virtual NumericalIE<T>* Clone(        ) const;

private:
    void Normalize(SMat& m) const
    {
        auto nn=OuterProduct(itsNormalizations);
        m=DirectMultiply(m,nn);
//        ::Normalize(m,itsNormalizations);
    }
    void Normalize(Vec& v) const
    {
        v=DirectMultiply(v,itsNormalizations);
    }
    void Normalize(const RVec&, Mat&, const RVec&) const;

    
    void    CheckInitialized       () const;
    const   RVec GetNumericalNormalization(bs_t&) const;
    virtual SMat MakeOverlap  (const TBasisFunction<T>& ) const;
    virtual SMat MakeRepulsion(const TBasisFunction<T>& ) const;

    IDRef<bs_t > itsBasisSet;       //Permenant hook into the basis set.
    RVec                       itsNormalizations; //Normalization constants.

    Mesh*              itsMesh;
    MeshIntegrator<T>* itsIntegrator;
};

#endif //_NumricalIEImp_H_
