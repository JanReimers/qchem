// File: HeapDB.H  Implement an integral data base that uses heap for storage.
#ifndef _HeapDB_H_
#define _HeapDB_H_



#include "IntegralDataBase.H"
#include "Misc/Persistent/IDRef.H"
#include "Misc/ERI4.H"
#include "Misc/rc_ptr.H"

//-----------------------------------------------------------------------------
//
//  This class implements an IntegralDataBase that simply stores all integrals
//  on the heap.  This means RAM and/or vitual memory in a real OS.
//
template <class T> class HeapDB :
    virtual public IntegralDataBase<T>
{
protected:
    typedef IntegralDataBase<T> Base;
    typedef typename Base::Vec   Vec;
    typedef typename Base::Mat   Mat;
    typedef typename Base::RVec  RVec;
    typedef typename Base::SMat  SMat;
    typedef typename Base::ERI3 ERI3;
public:

    HeapDB();
    virtual void  Insert(const TBasisSet<T>*,const NumericalIE<T>*);
    virtual void  Insert(const TBasisSet<T>*,const AnalyticIE<T>*);
    virtual void  Insert(const BasisGroup*);
    virtual void  Insert(const ERI4& J,const ERI4& K);
    virtual bool  operator==(const IntegralDataBase<T>&) const;

    virtual const SMat&    GetOverlap    (                             );
    virtual const Mat&     GetOverlap    (const TBasisSet<T>&          );
    virtual const Vec&     GetOverlap    (const ScalarFunction<double>&);
    virtual const ERI3&    GetOverlap3C  (const TBasisSet<T>&          );

    virtual const SMat&    GetRepulsion  (                             );
    virtual const Mat&     GetRepulsion  (const TBasisSet<T>&          );
    virtual const Vec&     GetRepulsion  (const ScalarFunction<double>&);
    virtual const ERI3&    GetRepulsion3C(const TBasisSet<T>&          );
    virtual       ERI4view GetRepulsion4C(const TBasisSet<T>* bs_cd    );
    virtual       ERI4view GetExchange4C (const TBasisSet<T>* bs_cd    );

    virtual const RVec&    GetNormalization (                );
    virtual const RVec&    GetCharge        (                );

    virtual const SMat&    GetKinetic         (              );
    virtual const SMat&    GetNuclear         (const Cluster&);
    virtual const SMat&    GetInverseOverlap  (              );
    virtual const SMat&    GetInverseRepulsion(              );

    virtual const NumericalIE<T>*  GetIntegralEngine(              ) const;

    virtual std::ostream&  Write(std::ostream&) const;
    virtual std::istream&  Read (std::istream&)      ;
    virtual IntegralDataBase<T>* Clone() const;

private:
    using StreamableObject::Binary;
    using StreamableObject::Ascii;
    using StreamableObject::Pretty;

    void WipeCleanAllData   ()      ;
    void BuildERIs();

    typedef std::vector<UniqueID::IDtype>      IDList;

    const NumericalIE<T>* itsNumericalIE;
    const AnalyticIE<T>*     itsAnalyticIE; //Integral engine for analytic integrals.
    const BasisSet*          itsBasisSet;
    const BasisGroup*        itsBasisGroup;

    SMat itsSelfOverlap  ; //Overlap with its self.
    Vec  itsRSFOverlap   ; //Overlap with ScalarFunction.
    SMat itsSelfRepulsion; //Repulsion with its self.
    Vec  itsRSFRepulsion ; //Repulsion with ScalarFunction.
    SMat itsSelfKinetic  ; //Self kinetic energy matrix.
    RVec itsNormalization; //Normalization constants.
    RVec itsCharge       ; //Integrated charge.
    SMat itsInvOverlap   ; //Inverse of the overlap matrix.
    SMat itsInvRepulsion ; //Inverse of the Repulsion matrix.

    bool itsSelfOverlapFlag     ; //These are set to true if they have been calculated.
    bool itsSelfRepulsionFlag   ;
    bool itsSelfKineticFlag     ;
    bool itsNormalizationFlag   ;
    bool itsChargeFlag          ;
    bool itsInvOverlapFlag      ;
    bool itsInvRepulsionFlag    ;

    Mat     itsBSOverlaps  ;       //Cross overlap with one other basis set.
    Mat     itsBSRepulsions;       //Cross repulsion with one other basis set.
    ERI3   itsNuclears    ;       //Nuclear attraction with a cluster of atoms.
    ERI3   its3CenterOverlaps;    //Overlap   with a whole basis set.  Only one basis is buffered.
    ERI3   its3CenterRepulsions;  //Repulsion with a whole basis set.  Only one basis is buffered.

    UniqueID::IDtype itsOverlapBasisSets;      //Basis set for which cross overlap has been calculated.
    UniqueID::IDtype itsRepulsionBasisSets;    //Basis set for which cross repulsion has been calculated.
    IDList           itsNuclearClusters;       //List of clusters for which nuclear attractions have been calculated.
    UniqueID::IDtype its3CenterOverlapBS;      //Basis set for which 3 center overlap has been calculated.
    UniqueID::IDtype its3CenterRepulsionBS;    //Basis set for which 3 center repulsion has been calculated.

    ERI4 itsJTable; //Coulomb integrals.
    ERI4 itsKTable; //Exchange integrals.
//    std::map<
};

#endif
