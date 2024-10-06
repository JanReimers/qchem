// File: dbheap.cpp  Implement a heap storage integral data base.



#include "DFTDataBase/HeapDB/HeapDB.H"
#include "QuantumNumber.H"
#include "BasisSet.H"
#include "Cluster.H"
#include "IntegralEngine1.H"
#include "Misc/ERIList.H"
#include "Misc/ERIProxy.H"
#include "oml/vector.h"
#include "oml/imp/binio.h"
#include "oml/vector.h"
#include "Misc/stl_io.h"
#include <iostream>
#include <cassert>
#include <algorithm> //find
#include <stdlib.h>

//------------------------------------------------------------------------
//
//  Construction zone.
//
template <class T> HeapDB<T>::HeapDB()
    :itsIntegralEngine    (0)
    ,itsIntegralEngine1   (0)
    ,itsBasisSet          (0)
    ,itsSelfOverlapFlag   (false)
    ,itsSelfRepulsionFlag (false)
    ,itsSelfKineticFlag   (false)
    ,itsNormalizationFlag (false)
    ,itsChargeFlag        (false)
    ,itsInvOverlapFlag    (false)
    ,itsInvRepulsionFlag  (false)
    ,itsBSOverlaps        ()
    ,itsBSRepulsions      ()
    ,itsNuclears          (new MList)
    ,its3CenterOverlaps   (new MList)
    ,its3CenterRepulsions (new MList)
    ,its4CenterRepulsions ()
    ,its4CenterExchange   ()
    ,itsOverlapBasisSets  (0)
    ,itsRepulsionBasisSets(0)
    ,its3CenterOverlapBS  (0)
    ,its3CenterRepulsionBS(0)
{}

// TODO Why do we need a copy constructor?
template <class T> HeapDB<T>::HeapDB(const HeapDB&)
    :itsIntegralEngine    (0)
    ,itsIntegralEngine1   (0)
    ,itsBasisSet          (0)
    ,itsSelfOverlapFlag   (false)
    ,itsSelfRepulsionFlag (false)
    ,itsSelfKineticFlag   (false)
    ,itsNormalizationFlag (false)
    ,itsChargeFlag        (false)
    ,itsInvOverlapFlag    (false)
    ,itsInvRepulsionFlag  (false)
    ,itsBSOverlaps        ()
    ,itsBSRepulsions      ()
    ,itsNuclears          (new MList)
    ,its3CenterOverlaps   (new MList)
    ,its3CenterRepulsions (new MList)
    ,its4CenterRepulsions ()
    ,its4CenterExchange   ()
    ,itsOverlapBasisSets  (0)
    ,itsRepulsionBasisSets(0)
    ,its3CenterOverlapBS  (0)
    ,its3CenterRepulsionBS(0)
{
     assert(false);
}

template <class T>HeapDB<T>& HeapDB<T>::operator=(const HeapDB<T>& other)
{
    assert(false);
    return *this;
}


template <class T> void HeapDB<T>::WipeCleanAllData()
{
    itsSelfOverlap  .SetLimits(MatLimits());
    itsRSFOverlap   .SetLimits(VecLimits());
    itsSelfRepulsion.SetLimits(MatLimits());
    itsRSFRepulsion .SetLimits(VecLimits());
    itsSelfKinetic  .SetLimits(MatLimits());
    itsNormalization.SetLimits(VecLimits());
    itsCharge       .SetLimits(VecLimits());
    itsInvOverlap   .SetLimits(MatLimits());
    itsInvRepulsion .SetLimits(MatLimits());
    itsBSOverlaps   .SetLimits(MatLimits());
    itsBSRepulsions .SetLimits(MatLimits());

    itsSelfOverlapFlag  =false;
    itsSelfRepulsionFlag=false;
    itsSelfKineticFlag  =false;
    itsNormalizationFlag=false;
    itsChargeFlag       =false;
    itsInvOverlapFlag   =false;
    itsInvRepulsionFlag =false;

    itsNuclears         ->Empty();
    its3CenterOverlaps  ->Empty();
    its3CenterRepulsions->Empty();
    its4CenterRepulsions .Empty();  //TODO We need to decide who owns the ERIList. There are many shallow copies.
    its4CenterExchange   .Empty();  //TODO We need to decide who owns the ERIList. There are many shallow copies.

    itsOverlapBasisSets  =0;
    itsRepulsionBasisSets=0;
    itsNuclearClusters   .clear();

    its3CenterOverlapBS  =0;
    its3CenterRepulsionBS=0;

}

template <class T> void HeapDB<T>::Insert(const TBasisSet<T>* bs,const IntegralEngine<T>* ie)
{
    assert(bs);
    assert(ie);
    itsBasisSet=bs;
    itsIntegralEngine=ie;
}
template <class T> void HeapDB<T>::Insert(const IntegralEngine1<T>* ie1)
{
    assert(ie1);
    itsIntegralEngine1=ie1;
}

template <class T> void HeapDB<T>::Insert(const BasisGroup* bg)
{
    assert(bg);
    itsBasisGroup=bg;
}

template <class T> void HeapDB<T>::Insert(const ERIList& Coulomb, const ERIList& exchange)
{
    its4CenterRepulsions=Coulomb;
    its4CenterExchange=exchange;
}
template <class T> void HeapDB<T>::Insert(const ERIList1& J, const ERIList1& K)
{
    itsJTable=J;
    itsKTable=K;
}

template <class T> bool HeapDB<T>::operator==(const IntegralDataBase<T>& idb) const
{
    if (!itsBasisSet || itsBasisSet->GetID()==0) return false;
    const HeapDB* hdb=dynamic_cast<const HeapDB*>(&idb);
    if (!hdb) return false;
    if (!hdb->itsBasisSet || hdb->itsBasisSet->GetID()==0) return false;
    return itsBasisSet->GetID()==hdb->itsBasisSet->GetID();
}
//-------------------------------------------------------------------------
//
//  Assorted fluff.
//
template <class T> const IntegralEngine<T>* HeapDB<T>::GetIntegralEngine() const
{
    assert(itsIntegralEngine);
    return itsIntegralEngine;
}

template <class T> IntegralDataBase<T>* HeapDB<T>::Clone() const
{
    return new HeapDB(*this);
}


//------------------------------------------------------------------
//
//  Pickle and un pickle the integral engine, all matricies, and
//  all lists of ID's.
//
template <class T> std::ostream& HeapDB<T>::Write(std::ostream& os) const
{
//
//  Write the IntegralEngine, which better be initialized.
//
    if(!Binary()) os << std::endl;
    UniqueID::Write(os);

//
//  Write the overlap, kinetic, normalization and charge if they are defined.
//
    os << itsSelfOverlap
    << itsRSFOverlap
    << itsSelfRepulsion
    << itsRSFRepulsion
    << itsSelfKinetic
    << itsNormalization
    << itsCharge
    << itsInvOverlap
    << itsInvRepulsion
    << itsBSOverlaps
    << itsBSRepulsions
    ;
//
//  Write all the big list of matricies.
//
    os << *itsNuclears
    << *its3CenterOverlaps
    << *its3CenterRepulsions
    //<< its4CenterRepulsions
    ;

    os << itsNuclearClusters;

    if (Binary())
    {
        BinaryWrite(itsSelfOverlapFlag   ,os);
        BinaryWrite(itsSelfRepulsionFlag ,os);
        BinaryWrite(itsSelfKineticFlag   ,os);
        BinaryWrite(itsNormalizationFlag ,os);
        BinaryWrite(itsChargeFlag        ,os);
        BinaryWrite(itsInvOverlapFlag    ,os);
        BinaryWrite(itsInvRepulsionFlag  ,os);
        BinaryWrite(itsOverlapBasisSets  ,os);
        BinaryWrite(itsRepulsionBasisSets,os);
        BinaryWrite(its3CenterOverlapBS  ,os);
        BinaryWrite(its3CenterRepulsionBS,os);

    }
    else
    {
        os << itsSelfOverlapFlag  << " " << itsSelfRepulsionFlag  << " "
        << itsSelfKineticFlag  << " " << itsNormalizationFlag  << " "
        << itsChargeFlag       << " " << itsInvOverlapFlag     << " "
        << itsInvRepulsionFlag << " ";
        os << itsOverlapBasisSets << " " << itsRepulsionBasisSets << " "
        << its3CenterOverlapBS << " " << its3CenterRepulsionBS << " ";
    }
    return os;
}

template <class T> std::istream& HeapDB<T>::Read (std::istream& is)
{
    UniqueID::Read(is);

//
//  For the next four matricies, the bool flag must be check to see
//  if the matrix is really defined.
//
    is >> itsSelfOverlap
    >> itsRSFOverlap
    >> itsSelfRepulsion
    >> itsRSFRepulsion
    >> itsSelfKinetic
    >> itsNormalization
    >> itsCharge
    >> itsInvOverlap
    >> itsInvRepulsion
    >> itsBSOverlaps
    >> itsBSRepulsions
    ;

//
//  Read all the big list of matricies.
//
    is >> *itsNuclears
    >> *its3CenterOverlaps
    >> *its3CenterRepulsions
    ;


    is >> itsNuclearClusters         ;

    if (Binary())
    {
        BinaryRead(itsSelfOverlapFlag   ,is);
        BinaryRead(itsSelfRepulsionFlag ,is);
        BinaryRead(itsSelfKineticFlag   ,is);
        BinaryRead(itsNormalizationFlag ,is);
        BinaryRead(itsChargeFlag        ,is);
        BinaryRead(itsInvOverlapFlag    ,is);
        BinaryRead(itsInvRepulsionFlag  ,is);
        BinaryRead(itsOverlapBasisSets  ,is);
        BinaryRead(itsRepulsionBasisSets,is);
        BinaryRead(its3CenterOverlapBS  ,is);
        BinaryRead(its3CenterRepulsionBS,is);

    }
    else
    {
        is
        >> itsSelfOverlapFlag  >> itsSelfRepulsionFlag
        >> itsSelfKineticFlag  >> itsNormalizationFlag
        >> itsChargeFlag       >> itsInvOverlapFlag
        >> itsInvRepulsionFlag;
        is >> itsOverlapBasisSets >> itsRepulsionBasisSets >> its3CenterOverlapBS >> its3CenterRepulsionBS;
        is.get();
    }
    return is;
}





//---------------------------------------------------------------------------------
//
//  Calculate some matricies for which no lists need to be maintained.
//
template <class T>  const typename HeapDB<T>::SMat& HeapDB<T>::GetOverlap()
{
    assert(itsIntegralEngine);
    if(!itsSelfOverlapFlag)
    {
        if (itsIntegralEngine1)
            itsSelfOverlap=itsIntegralEngine1->MakeOverlap();
        else
            itsSelfOverlap=itsIntegralEngine->MakeOverlap();
        itsSelfOverlapFlag=true;
    }
    return itsSelfOverlap;
}

template <class T> const typename HeapDB<T>::Mat& HeapDB<T>::GetOverlap(const TBasisSet<T>& theBasisSet)
{
    // No UT coverage.
    assert(itsIntegralEngine);
    if (itsOverlapBasisSets!=theBasisSet.GetID())
    {
//        if (itsIntegralEngine1)
//            itsBSOverlaps = itsIntegralEngine1->MakeOverlap(theBasisSet.GetIntegralEngine1());
//        else
            itsBSOverlaps = itsIntegralEngine->MakeOverlap(theBasisSet);
        itsOverlapBasisSets = theBasisSet.GetID();
    }
    return itsBSOverlaps;
}

template <class T> const typename HeapDB<T>::Vec& HeapDB<T>::GetOverlap(const ScalarFunction<double>& f)
{
    //Only used for numerical IE.
    assert(itsIntegralEngine);
    itsRSFOverlap=itsIntegralEngine->MakeOverlap(f);

    return itsRSFOverlap;
}



template <class T> const typename HeapDB<T>::SMat& HeapDB<T>::GetRepulsion()
{
    assert(itsIntegralEngine);
    if(!itsSelfRepulsionFlag)
    {
        if (itsIntegralEngine1)
            itsSelfRepulsion=itsIntegralEngine1->MakeRepulsion();
        else
            itsSelfRepulsion=itsIntegralEngine->MakeRepulsion();
        itsSelfRepulsionFlag=true;
    }
    return itsSelfRepulsion;
}

template <class T> const typename HeapDB<T>::Vec& HeapDB<T>::GetRepulsion(const ScalarFunction<double>& f)
{
    assert(itsIntegralEngine);
    itsRSFRepulsion=itsIntegralEngine->MakeRepulsion(f);
    return itsRSFRepulsion;
}

template <class T> const typename HeapDB<T>::Mat& HeapDB<T>::GetRepulsion(const TBasisSet<T>& theBasisSet)
{
    assert(itsIntegralEngine);
    if (itsRepulsionBasisSets!=theBasisSet.GetID())
    {
//        if (itsIntegralEngine1)
//            itsBSRepulsions = itsIntegralEngine1->MakeRepulsion(theBasisSet.GetIntegralEngine1());
//        else
            itsBSRepulsions = itsIntegralEngine->MakeRepulsion(theBasisSet);
        itsRepulsionBasisSets = theBasisSet.GetID();
    }
    return itsBSRepulsions;
}

template <class T> const typename HeapDB<T>::SMat& HeapDB<T>::GetKinetic()
{
    assert(itsIntegralEngine);
    if(!itsSelfKineticFlag)
    {
         if (itsIntegralEngine1)
            itsSelfKinetic=itsIntegralEngine1->MakeKinetic();
        else
            itsSelfKinetic=itsIntegralEngine->MakeKinetic();
        itsSelfKineticFlag=true;
    }
    return itsSelfKinetic;
}

template <class T> const typename HeapDB<T>::RVec& HeapDB<T>::GetNormalization()
{
    assert(itsIntegralEngine);
    if(!itsNormalizationFlag)
    {
//        if (itsIntegralEngine1)
//            itsNormalization=itsIntegralEngine1->MakeNormalization();
//        else
            itsNormalization=itsIntegralEngine->MakeNormalization();
        itsNormalizationFlag=true;
    }
    return itsNormalization;
}

template <class T> const typename HeapDB<T>::RVec& HeapDB<T>::GetCharge()
{
    assert(itsIntegralEngine);
    if(!itsChargeFlag)
    {
//        if (itsIntegralEngine1)
//            itsCharge=itsIntegralEngine1->MakeCharge();
//        else
            itsCharge=itsIntegralEngine->MakeCharge();
        itsChargeFlag=true;
    }
    return itsCharge;
}

template <class T> const typename HeapDB<T>::MList& HeapDB<T>::GetOverlap3C(const TBasisSet<T>& bs)
{
    assert(itsIntegralEngine);
    if (its3CenterOverlapBS!=bs.GetID())
    {
        if (itsIntegralEngine1)
            itsIntegralEngine1->MakeOverlap3C(*its3CenterOverlaps,bs.GetIntegralEngine1());
        else
            itsIntegralEngine->MakeOverlap3C(*its3CenterOverlaps,bs);
        its3CenterOverlapBS=bs.GetID();
    }
    return *its3CenterOverlaps;
}

template <class T> const typename HeapDB<T>::MList& HeapDB<T>::GetRepulsion3C(const TBasisSet<T>& bs)
{
    assert(itsIntegralEngine);
    if (its3CenterRepulsionBS!=bs.GetID())
    {
//        if (itsIntegralEngine1)
//            itsIntegralEngine1->MakeRepulsion3C(*its3CenterRepulsions,bs.GetIntegralEngine1());
//        else
            itsIntegralEngine->MakeRepulsion3C(*its3CenterRepulsions,bs);
        its3CenterRepulsionBS=bs.GetID();
    }
    return *its3CenterRepulsions;
}


template <class T> const ERIProxy HeapDB<T>::GetRepulsion4C(const TBasisSet<T>* bs_cd)
{
    assert(itsBasisSet);
    assert(bs_cd);
//    std::cout << "HeapDB<T>::GetRepulsion4C" << std::endl;
//    std::cout << "Repulsions: " << itsBasisSet->GetID() << " " << bs_cd->GetID()  << std::endl;
//    std::cout << "Repulsions: " << itsBasisSet->GetQuantumNumber() << " " << bs_cd->GetQuantumNumber()  << std::endl;

    if (its4CenterRepulsions.GetSize()==0)
    {
        assert(itsBasisGroup);
//        if (itsIntegralEngine1)
//            itsIntegralEngine1->MakeRepulsion4C(its4CenterRepulsions,its4CenterExchange,itsBasisGroup->Flatten());
//        else
            itsIntegralEngine->MakeRepulsion4C(its4CenterRepulsions,its4CenterExchange,itsBasisGroup);
        //
        //  This is how we pass the ERI tables to all the other DBs.
        //
        itsBasisGroup->Insert(its4CenterRepulsions,its4CenterExchange); //Eeach basis set has its own HeapDB.
     
    }
    
    return ERIProxy(its4CenterRepulsions,itsBasisSet->GetStartIndex(),bs_cd->GetStartIndex());
}

template <class T> const ERIProxy HeapDB<T>::GetExchange4C(const TBasisSet<T>* bs_cd)
{
    assert(itsBasisSet);
    assert(bs_cd);
    //std::cout << "Exchange: " << itsBasisSet->GetQuantumNumber() << " " << bs_cd->GetQuantumNumber()  << std::endl;

    if (its4CenterExchange.GetSize()==0)
        return ERIProxy(its4CenterRepulsions,itsBasisSet->GetStartIndex(),bs_cd->GetStartIndex());
    else
        return ERIProxy(its4CenterExchange,itsBasisSet->GetStartIndex(),bs_cd->GetStartIndex());

}

template <class T> void HeapDB<T>::BuildERIs()
{
    assert(itsIntegralEngine1);
    assert(itsJTable.GetSize()==0);
    assert(itsKTable.GetSize()==0); //They should be synchronized.
    auto [J,K]=itsIntegralEngine1->Make4C(itsBasisGroup->Flatten());
    itsBasisGroup->Insert(J,K); //Eeach basis set has its own HeapDB.
//    index_t N=bs_cd->GetNumFunctions();
//    for (int ia=1;ia<=N;ia++)
//    for (int ib=1;ib<=N;ib++)
//    for (int ic=1;ic<=N;ic++)
//    for (int id=1;id<=N;id++)
//    {
//    //                std::cout << "J(" << ia << "," << ib << "," << ic << "," << id << ")   " 
//    //                << J.GetIndex(ia,ib,ic,id) << " " << J2.GetIndex(ia,ib,ic,id) << " " << J(ia,ib,ic,id) << " " << J2(ia,ib,ic,id) << std::endl;
//               assert(itsJTable(ia,ib,ic,id)==its4CenterRepulsions(ia,ib,ic,id));
//    //                std::cout << "K(" << ia << "," << ib << "," << ic << "," << id << ")   " 
//    //                << K.GetIndex(ia,ib,ic,id) << " " << K2.GetIndex(ia,ib,ic,id) << " " << K(ia,ib,ic,id) << " " << K2(ia,ib,ic,id) << std::endl;
//                assert(itsKTable(ia,ib,ic,id)==its4CenterExchange(ia,ib,ic,id));
//            }
//        }
   
}
template <class T> ERIProxy1 HeapDB<T>::GetRepulsion4C_1(const TBasisSet<T>* bs_cd)
{
   if (itsJTable.GetSize()==0) BuildERIs(); 
   return ERIProxy1(itsJTable,itsBasisSet->GetStartIndex(),bs_cd->GetStartIndex());
}

template <class T> ERIProxy1 HeapDB<T>::GetExchange4C_1 (const TBasisSet<T>* bs_cd)
{
   if (itsJTable.GetSize()==0) BuildERIs(); 
   return ERIProxy1(itsKTable,itsBasisSet->GetStartIndex(),bs_cd->GetStartIndex());
}
//-------------------------------------------------------------------------
//
//  These guys have to check the ID lists to see if it integrals are
//  already known.
//

template <class T> const typename HeapDB<T>::SMat& HeapDB<T>::GetNuclear(const Cluster& theCluster)
{
    assert(itsIntegralEngine);
//    std::cout << theCluster.GetID() << std::endl;
//    std::cout << "itsNuclearClusters.size()=" << itsNuclearClusters.size() << std::endl;
    unsigned int index=itsNuclearClusters.size();
    if (index>0)
    {
        auto i=std::find(itsNuclearClusters.begin(),itsNuclearClusters.end(),theCluster.GetID());
        if (i!=itsNuclearClusters.end()) index=i-itsNuclearClusters.begin();
    }
//    std::cout << "HeapDB<T>::GetNuclear index=" << index <<  std::endl;
    if (index == itsNuclearClusters.size())
    {
//        std::cout << "Creating cluster index=" << index << " ID=" << theCluster.GetID() << std::endl;
        if (itsIntegralEngine1)
            itsNuclears->Add(itsIntegralEngine1->MakeNuclear(theCluster));
        else
            itsNuclears->Add(itsIntegralEngine->MakeNuclear(theCluster));
        itsNuclearClusters.push_back(theCluster.GetID());
    }
    return (*itsNuclears)[index];
}


template <class T> const typename HeapDB<T>::SMat& HeapDB<T>::GetInverseOverlap()
{
    assert(itsIntegralEngine);
    if (!itsInvOverlapFlag)
    {
        itsInvOverlap = itsIntegralEngine->MakeInverse(GetOverlap());
        itsInvOverlapFlag = true;
    }
    return itsInvOverlap;
}

template <class T> const typename HeapDB<T>::SMat& HeapDB<T>::GetInverseRepulsion()
{
    assert(itsIntegralEngine);
    if (!itsInvRepulsionFlag)
    {
        SMat repulsion=GetRepulsion();
//        std::cout << "repulsion=" << repulsion << std::endl;
        itsInvRepulsion = itsIntegralEngine->MakeInverse(repulsion);
//        std::cout << "itsInvRepulsion=" << itsInvRepulsion << std::endl;
        itsInvRepulsionFlag = true;
    }
    return itsInvRepulsion;
}

template <> const HeapDB<std::complex<double> >::SMat& HeapDB<std::complex<double> >::
GetInverseOverlap()
{
    std::cerr << "Sorry, inverse of complex matrix is not implemented" << std::endl;
    exit(-1);
    return itsInvOverlap;
}

template <> const HeapDB<std::complex<double> >::SMat& HeapDB<std::complex<double> >::
GetInverseRepulsion()
{
    std::cerr << "Sorry, inverse of complex matrix is not implemented" << std::endl;
    exit(-1);
    return itsInvRepulsion;
}

template class HeapDB<double>;
template class HeapDB<std::complex<double> >;
