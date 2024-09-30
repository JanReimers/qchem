// File: HeapDB.H  Implement an integral data base that uses heap for storage.
#ifndef _HeapDB_H_
#define _HeapDB_H_



#include "BasisSet/IntegralEngine.H"
#include "BasisSet/IntegralDataBase.H"
#include "Misc/Persistent/IDRef.H"
#include "Misc/ERIList.H"
#include "oml/vector.h"
#include "oml/matrix.h"
#include "oml/smatrix.h"
#include "Misc/MatrixList.H"
#include "Misc/rc_ptr.H"


class BasisSet;

//-----------------------------------------------------------------------------
//
//  This class implements an IntegralDataBase that simply stores all integrals
//  on the heap.  This means RAM and/or vitual memory in a real OS.
//
template <class T> class HeapDB :
    virtual public IntegralDataBase<T>
{
protected:
    typedef IntegralEngine<T> Base;
    typedef typename Base::Vec   Vec;
    typedef typename Base::Mat   Mat;
    typedef typename Base::RVec  RVec;
    typedef typename Base::SMat  SMat;
    typedef typename Base::MList MList;
public:
    HeapDB();
    HeapDB(const HeapDB&);
    HeapDB& operator=(const HeapDB&);
    virtual void  Insert(const TBasisSet<T>*,  const IntegralEngine<T>*);
    virtual void  Insert(const BasisGroup*);
    virtual void  Insert(const ERIList& c_eris,const ERIList& x_eris);
    virtual bool  operator==(const IntegralDataBase<T>&) const;

    virtual const SMat&    GetOverlap    (                             );
    virtual const Mat&     GetOverlap    (const TBasisSet<T>&          );
    virtual const Vec&     GetOverlap    (const ScalarFunction<double>&);
    virtual const MList&   GetOverlap3C  (const TBasisSet<T>&          );

    virtual const SMat&    GetRepulsion  (                             );
    virtual const Mat&     GetRepulsion  (const TBasisSet<T>&          );
    virtual const Vec&     GetRepulsion  (const ScalarFunction<double>&);
    virtual const MList&   GetRepulsion3C(const TBasisSet<T>&          );
    virtual const ERIProxy GetRepulsion4C(const TBasisSet<T>* bs_cd    );
    virtual const ERIProxy GetExchange4C (const TBasisSet<T>* bs_cd    );

    virtual const RVec&    GetNormalization (                );
    virtual const RVec&    GetCharge        (                );

    virtual const SMat&    GetKinetic         (              );
    virtual const SMat&    GetNuclear         (const Cluster&);
    virtual const SMat&    GetInverseOverlap  (              );
    virtual const SMat&    GetInverseRepulsion(              );

    virtual const IntegralEngine<T>*  GetIntegralEngine(              ) const;

    virtual std::ostream&  Write(std::ostream&) const;
    virtual std::istream&  Read (std::istream&)      ;
    virtual IntegralDataBase<T>* Clone() const;

private:
    using StreamableObject::Binary;
    using StreamableObject::Ascii;
    using StreamableObject::Pretty;

    void WipeCleanAllData   ()      ;

    typedef std::vector<UniqueID::IDtype>      IDList;

    const IntegralEngine<T>* itsIntegralEngine;
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
    rc_ptr<MList>   itsNuclears    ;       //Nuclear attraction with a cluster of atoms.
    rc_ptr<MList>   its3CenterOverlaps;    //Overlap   with a whole basis set.  Only one basis is buffered.
    rc_ptr<MList>   its3CenterRepulsions;  //Repulsion with a whole basis set.  Only one basis is buffered.
    ERIList   its4CenterRepulsions; //Supermatrix of 4C repulsion integrals.
    ERIList   its4CenterExchange; //Supermatrix of 4C exchange repulsion integrals.  Only used for certain basis sets.

    UniqueID::IDtype itsOverlapBasisSets;      //Basis set for which cross overlap has been calculated.
    UniqueID::IDtype itsRepulsionBasisSets;    //Basis set for which cross repulsion has been calculated.
    IDList           itsNuclearClusters;       //List of clusters for which nuclear attractions have been calculated.
    UniqueID::IDtype its3CenterOverlapBS;      //Basis set for which 3 center overlap has been calculated.
    UniqueID::IDtype its3CenterRepulsionBS;    //Basis set for which 3 center repulsion has been calculated.
};

#endif
