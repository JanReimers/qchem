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

public:
    NumericalIEImp(                  );
    NumericalIEImp(const NumericalIEImp&);
    NumericalIEImp(Mesh*             );
    ~NumericalIEImp(                  );

    virtual void Insert(const TBasisSet<T>*);

    virtual SMat MakeOverlap    (                             ) const;
    virtual Mat  MakeOverlap    (const TBasisSet<T>&          ) const;
    virtual Vec  MakeOverlap    (const ScalarFunction<double>&) const;

    virtual SMat MakeRepulsion  (                             ) const;
    virtual Mat  MakeRepulsion  (const TBasisSet<T>&          ) const;
    virtual Vec  MakeRepulsion  (const ScalarFunction<double>&) const;

    virtual RVec MakeNormalization(                 ) const;
    virtual RVec MakeCharge       (                 ) const;
#ifdef USE_FOR_DEBUGGING_ANALYTIC
    virtual SMat MakeKinetic  (                     ) const;
    virtual SMat MakeNuclear  (const Cluster&       ) const;
    virtual void MakeOverlap3C  (MList&,const TBasisSet<T>&   ) const;
    virtual void MakeRepulsion3C(MList&,const TBasisSet<T>&   ) const;
    virtual void MakeRepulsion4C(ERIList&,ERIList&,const BasisGroup*   ) const;
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
    const   RVec GetNumericalNormalization(const TBasisSet<T>&) const;
    virtual SMat MakeOverlap  (const TBasisFunction<T>& ) const;
    virtual SMat MakeRepulsion(const TBasisFunction<T>& ) const;

    IDRef<const TBasisSet<T> > itsBasisSet;       //Permenant hook into the basis set.
    RVec                       itsNormalizations; //Normalization constants.

    Mesh*              itsMesh;
    MeshIntegrator<T>* itsIntegrator;
};

#endif //_NumricalIEImp_H_
