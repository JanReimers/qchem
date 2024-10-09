// File: NumricalIEImp.H  Numerical integral engines.
#ifndef _NumricalIEImp_H_
#define _NumricalIEImp_H_



#include "NumericalIE.H"

template <class T> class MeshIntegrator;

template <class T> class NumericalIEImp
    : public virtual NumericalIE<T>
{
    typedef NumericalIE<T> Base;
    typedef typename Base::Vec   Vec;
    typedef typename Base::Mat   Mat;
    typedef typename Base::RVec  RVec;
    typedef typename Base::SMat  SMat;
    typedef typename Base::bs_t bs_t;
    typedef typename Base::Rf Rf;
    typedef typename Base::Vf Vf;
public:
    NumericalIEImp(                     );
    NumericalIEImp(const NumericalIEImp&);
    NumericalIEImp(Mesh*                );
    ~NumericalIEImp(                    );

    virtual RVec MakeNormalization(bs_t& a) const;
    virtual Vec  MakeOverlap    (Vf&    , Rf&) const;
    virtual Vec  MakeRepulsion  (Vf&    , Rf&) const;
    //
    //  These currently have no UT coverage.
    //
#ifndef UT_COVERAGE_ONLY
    virtual RVec MakeCharge       (bs_t& a) const;
    virtual SMat MakeOverlap      (bs_t& a) const;
    virtual SMat MakeRepulsion    (bs_t& a) const;
    virtual Mat  MakeOverlap      (bs_t& a,bs_t& b) const;
    virtual Mat  MakeRepulsion    (bs_t& a,bs_t& b) const;
#endif



#ifdef USE_FOR_DEBUGGING_ANALYTIC
    virtual SMat MakeKinetic  (                     ) const;
    virtual SMat MakeNuclear  (const Cluster&       ) const;
    virtual void MakeOverlap3C  (ERI3&,bs_t&   ) const;
    virtual void MakeRepulsion3C(ERI3&,bs_t&   ) const;
    virtual void MakeRepulsion4C(ERI4&,ERIList&,const BasisGroup*   ) const;
#endif

    virtual NumericalIE<T>* Clone(        ) const;

private:
    const   RVec GetNumericalNormalization(bs_t&) const;
    
    Mesh*              itsMesh;
    MeshIntegrator<T>* itsIntegrator;
};

#endif //_NumricalIEImp_H_
