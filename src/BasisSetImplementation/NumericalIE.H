// File: NumricalIE.H  Numerical integral engines.
#ifndef _NumricalIE_H_
#define _NumricalIE_H_



#include "BasisSetImplementation/IntegralEngineImplementation.H"
#include "Mesh/MeshIntegrator.H"

template <class T> class NumericalIE
    : public virtual IntegralEngine<T>
    , private IntegralEngineImplementation<T>
{
    typedef IntegralEngineImplementation<T> Imp;
    typedef IntegralEngine<T> Base;
    typedef typename Base::Vec   Vec;
    typedef typename Base::Mat   Mat;
    typedef typename Base::RVec  RVec;
    typedef typename Base::SMat  SMat;
    typedef typename Base::MList MList;

public:
    NumericalIE(                  );
    NumericalIE(const NumericalIE&);
    NumericalIE(Mesh*             );
    ~NumericalIE(                  );

    virtual void Insert(const TBasisSet<T>*);

    virtual SMat MakeOverlap    (                             ) const;
    virtual Mat  MakeOverlap    (const TBasisSet<T>&          ) const;
    virtual Vec  MakeOverlap    (const ScalarFunction<double>&) const;
    virtual void MakeOverlap3C  (MList&,const TBasisSet<T>&   ) const;

    virtual SMat MakeRepulsion  (                             ) const;
    virtual Mat  MakeRepulsion  (const TBasisSet<T>&          ) const;
    virtual Vec  MakeRepulsion  (const ScalarFunction<double>&) const;
    virtual void MakeRepulsion3C(MList&,const TBasisSet<T>&   ) const;
    virtual void MakeRepulsion4C(ERIList&,ERIList&,const BasisGroup*   ) const;

    virtual RVec MakeNormalization(                 ) const;
    virtual RVec MakeCharge       (                 ) const;

    virtual SMat MakeKinetic  (                     ) const;
    virtual SMat MakeNuclear  (const Cluster&       ) const;

    virtual std::ostream&           Write(std::ostream&) const;
    virtual std::istream&           Read (std::istream&)      ;
    virtual IntegralEngine<T>* Clone(        ) const;

private:
    using Imp::CheckBasisSet;
    using Imp::Normalize;
    using Imp::itsNormalizations;
    using Imp::itsBasisSet;

    void    CheckInitialized       () const;
    const   RVec GetNumericalNormalization(const TBasisSet<T>&) const;
    virtual SMat MakeOverlap  (const TBasisFunction<T>& ) const;
    virtual SMat MakeRepulsion(const TBasisFunction<T>& ) const;

    Mesh*              itsMesh;
    MeshIntegrator<T>* itsIntegrator;
};

#endif //_NumricalIE_H_
