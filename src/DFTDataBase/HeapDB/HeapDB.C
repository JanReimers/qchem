// File: dbheap.cpp  Implement a heap storage integral data base.



#include "DFTDataBase/HeapDB/HeapDB.H"
//#include "QuantumNumber.H"
#include "BasisSet.H"
#include "Cluster.H"
#include "NumericalIE.H"
#include "AnalyticIE.H"
//#include "Misc/ERIList.H"
//#include "Misc/ERIProxy.H"
#include "oml/vector.h"
#include "Misc/stl_io.h"
#include <iostream>
//#include <cassert>
#include <algorithm> //find

//------------------------------------------------------------------------
//
//  Construction zone.
//
template <class T> HeapDB<T>::HeapDB()
    :itsNumericalIE    (0)
    ,itsAnalyticIE   (0)
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
    ,itsNuclears          ()
    ,its3CenterOverlaps   ()
    ,its3CenterRepulsions ()
    ,itsOverlapBasisSets  (0)
    ,itsRepulsionBasisSets(0)
    ,its3CenterOverlapBS  (0)
    ,its3CenterRepulsionBS(0)
{}



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

    itsNuclears         .clear();
    its3CenterOverlaps  .clear();
    its3CenterRepulsions.clear();

    itsOverlapBasisSets  =0;
    itsRepulsionBasisSets=0;
    itsNuclearClusters   .clear();

    its3CenterOverlapBS  =0;
    its3CenterRepulsionBS=0;

}

template <class T> void HeapDB<T>::Insert(const TBasisSet<T>* bs,const NumericalIE<T>* ie)
{
    assert(bs);
    assert(ie);
    itsBasisSet=bs;
    itsNumericalIE=ie;
}

template <class T> void HeapDB<T>::Insert(const TBasisSet<T>* bs,const AnalyticIE<T>* ie)
{
    assert(bs);
    assert(ie);
    itsBasisSet=bs;
    itsAnalyticIE=ie;
}

template <class T> void HeapDB<T>::Insert(const BasisGroup* bg)
{
    assert(bg);
    itsBasisGroup=bg;
}

//template <class T> void HeapDB<T>::Insert(const ERIList& Coulomb, const ERIList& exchange)
//{
//    its4CenterRepulsions=Coulomb;
//    its4CenterExchange=exchange;
//}


template <class T> void HeapDB<T>::Insert(const ERI4& J, const ERI4& K)
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
template <class T> const NumericalIE<T>* HeapDB<T>::GetIntegralEngine() const
{
    assert(itsNumericalIE);
    return itsNumericalIE;
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
    os << itsNuclears
    << its3CenterOverlaps
    << its3CenterRepulsions
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
    is >> itsNuclears
    >> its3CenterOverlaps
    >> its3CenterRepulsions
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
    assert(itsAnalyticIE || itsNumericalIE);
    if(!itsSelfOverlapFlag)
    {
        if (itsNumericalIE)
            itsSelfOverlap=itsNumericalIE->MakeOverlap();
        else
            itsSelfOverlap=itsAnalyticIE->MakeOverlap();
        itsSelfOverlapFlag=true;
    }
    return itsSelfOverlap;
}

template <class T> const typename HeapDB<T>::Mat& HeapDB<T>::GetOverlap(const TBasisSet<T>& theBasisSet)
{
    // No UT coverage.
    assert(false);
    assert(itsNumericalIE);
    if (itsOverlapBasisSets!=theBasisSet.GetID())
    {
        itsBSOverlaps = itsNumericalIE->MakeOverlap(theBasisSet);
        itsOverlapBasisSets = theBasisSet.GetID();
    }
    return itsBSOverlaps;
}

template <class T> const typename HeapDB<T>::Vec& HeapDB<T>::GetOverlap(const ScalarFunction<double>& f)
{
    //Only used for numerical IE.
    assert(itsNumericalIE);
    itsRSFOverlap=itsNumericalIE->MakeOverlap(f);

    return itsRSFOverlap;
}



template <class T> const typename HeapDB<T>::SMat& HeapDB<T>::GetRepulsion()
{
    assert(itsAnalyticIE || itsNumericalIE);
    if(!itsSelfRepulsionFlag)
    {
        if (itsNumericalIE)
            itsSelfRepulsion=itsNumericalIE->MakeRepulsion();
        else
            itsSelfRepulsion=itsAnalyticIE->MakeRepulsion();
        itsSelfRepulsionFlag=true;
    }
    return itsSelfRepulsion;
}

template <class T> const typename HeapDB<T>::Vec& HeapDB<T>::GetRepulsion(const ScalarFunction<double>& f)
{
    //Only used for numerical IE.
    assert(itsNumericalIE);
    itsRSFRepulsion=itsNumericalIE->MakeRepulsion(f);
    return itsRSFRepulsion;
}

template <class T> const typename HeapDB<T>::Mat& HeapDB<T>::GetRepulsion(const TBasisSet<T>& theBasisSet)
{
    assert(itsAnalyticIE || itsNumericalIE);
    if (itsRepulsionBasisSets!=theBasisSet.GetID())
    {
        if (itsNumericalIE)
            itsBSRepulsions = itsNumericalIE->MakeRepulsion(theBasisSet);
        else
            itsBSRepulsions = itsAnalyticIE->MakeRepulsion(theBasisSet.GetAnalyticIE());
        
        itsRepulsionBasisSets = theBasisSet.GetID();
    }
    return itsBSRepulsions;
}

template <class T> const typename HeapDB<T>::SMat& HeapDB<T>::GetKinetic()
{
    assert(itsAnalyticIE);
    if(!itsSelfKineticFlag)
    {
        itsSelfKinetic=itsAnalyticIE->MakeKinetic();
        itsSelfKineticFlag=true;
    }
    return itsSelfKinetic;
}

template <class T> const typename HeapDB<T>::RVec& HeapDB<T>::GetNormalization()
{
    assert(itsNumericalIE);
    if(!itsNormalizationFlag)
    {
        itsNormalization=itsNumericalIE->MakeNormalization();
        itsNormalizationFlag=true;
    }
    return itsNormalization;
}

template <class T> const typename HeapDB<T>::RVec& HeapDB<T>::GetCharge()
{
    assert(itsAnalyticIE || itsNumericalIE);
    if(!itsChargeFlag)
    {
        if (itsNumericalIE)
            itsCharge=itsNumericalIE->MakeCharge();
        else
            itsCharge=itsAnalyticIE->MakeCharge();

        itsChargeFlag=true;
    }
    return itsCharge;
}

template <class T> const typename HeapDB<T>::ERI3& HeapDB<T>::GetOverlap3C(const TBasisSet<T>& bs)
{
    assert(itsAnalyticIE);
    if (its3CenterOverlapBS!=bs.GetID())
    {
        itsAnalyticIE->MakeOverlap3C(its3CenterOverlaps,bs.GetAnalyticIE());
        its3CenterOverlapBS=bs.GetID();
    }
    return its3CenterOverlaps;
}

template <class T> const typename HeapDB<T>::ERI3& HeapDB<T>::GetRepulsion3C(const TBasisSet<T>& bs)
{
    assert(itsAnalyticIE);
    if (its3CenterRepulsionBS!=bs.GetID())
    {
        itsAnalyticIE->MakeRepulsion3C(its3CenterRepulsions,bs.GetAnalyticIE());
        its3CenterRepulsionBS=bs.GetID();
    }
    return its3CenterRepulsions;
}

template <class T> void HeapDB<T>::BuildERIs()
{
    assert(itsAnalyticIE);
    assert(itsJTable.GetSize()==0);
    assert(itsKTable.GetSize()==0); //They should be synchronized.
    auto [J,K]=itsAnalyticIE->Make4C(itsBasisGroup->Flatten());
    itsBasisGroup->Insert(J,K); //Eeach basis set has its own HeapDB.
}

template <class T> ERI4view HeapDB<T>::GetRepulsion4C(const TBasisSet<T>* bs_cd)
{
    assert(itsBasisSet);
    assert(itsAnalyticIE);
   if (itsJTable.GetSize()==0) BuildERIs(); 
   assert(bs_cd);
   int ss=bs_cd->GetStartIndex();
   return ERI4view(itsJTable,itsBasisSet->GetStartIndex(),ss);
}

template <class T> ERI4view HeapDB<T>::GetExchange4C (const TBasisSet<T>* bs_cd)
{
   if (itsJTable.GetSize()==0) BuildERIs(); 
   if (itsKTable.GetSize()==0)
        return ERI4view(itsJTable,itsBasisSet->GetStartIndex(),bs_cd->GetStartIndex());
    else
        return ERI4view(itsKTable,itsBasisSet->GetStartIndex(),bs_cd->GetStartIndex());
}
//-------------------------------------------------------------------------
//
//  These guys have to check the ID lists to see if it integrals are
//  already known.
//

template <class T> const typename HeapDB<T>::SMat& HeapDB<T>::GetNuclear(const Cluster& theCluster)
{
    assert(itsAnalyticIE);
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
        itsNuclears.push_back(itsAnalyticIE->MakeNuclear(theCluster));
        itsNuclearClusters.push_back(theCluster.GetID());
    }
    return itsNuclears[index];
}


template <class T> const typename HeapDB<T>::SMat& HeapDB<T>::GetInverseOverlap()
{
     //Only used for numerical IE.
    assert(itsNumericalIE);
    if (!itsInvOverlapFlag)
    {
        itsInvOverlap = itsNumericalIE->MakeInverse(GetOverlap());
        itsInvOverlapFlag = true;
    }
    return itsInvOverlap;
}

template <class T> const typename HeapDB<T>::SMat& HeapDB<T>::GetInverseRepulsion()
{
    assert(itsAnalyticIE);
    if (!itsInvRepulsionFlag)
    {
        SMat repulsion=GetRepulsion();
//        std::cout << "repulsion=" << repulsion << std::endl;
        itsInvRepulsion = itsAnalyticIE->MakeInverse(repulsion);
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
