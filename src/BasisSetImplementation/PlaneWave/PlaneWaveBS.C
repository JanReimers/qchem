// File: PlaneWaveBS.C  Polarized Gaussian basis set, for MO calculations.

#pragma implementation

#include "BasisSetImplementation/PlaneWave/PlaneWaveBF.H"
#include "BasisSetImplementation/PlaneWave/PlaneWaveBS.H"
#include "BasisSetImplementation/PlaneWave/PlaneWaveIE.H"
#include "BasisSetImplementation/PlaneWave/BlochQN.H"
#include "Cluster/Lattice.H"
#include "BasisSetImplementation/NumericalIE.H"
#include "BasisSet/IntegralDataBase.H"
#include "BasisSet/BasisSetBrowser.H"
#include "Cluster/Lattice.H"
#include "Misc/DFTDefines.H"
#include "oml/list.h"
#include <cassert>
#include <typeinfo>
#include <iostream>
#include <stdlib.h>

//#######################################################################
//
//  Concrete  gaussian basis set.
//
PlaneWaveBS::PlaneWaveBS()
    : BasisSetImplementation()
    , TBasisSetImplementation<std::complex<double> >()
{};

PlaneWaveBS::
PlaneWaveBS(IntegralDataBase<std::complex<double> >* theDB, const Cluster* theCluster, double Gmax, Mesh * theMesh)
    : BasisSetImplementation(new BlochQN(RVec3(0,0,0)))
    , TBasisSetImplementation<std::complex<double> >(theDB)
    , itsRLCell()
{
    const Lattice* DirLattice=dynamic_cast<const Lattice*>(theCluster);
    if(!DirLattice)
    {
        std::cerr << "Constructing plane waves requires a lattice object" << std::endl;
        std::cerr << "The cluster you supplied is of type " << typeid(*theCluster).name() << " not type Lattice" << std::endl;
        exit(-1);
    }


    double LatticeVolume=DirLattice->GetLatticeVolume();
    UnitCell DirCell=DirLattice->GetUnitCell();
    itsRLCell=DirCell.MakeReciprocalCell();

    List<RVec3> Basis;
    Basis.push_back(RVec3(0,0,0));
    Lattice RecLattice(itsRLCell,itsRLCell.GetNumCells(Gmax));
    RecLattice.Insert(new Atom);

    List<RVec3> itsGs= RecLattice.GetBondsInSphere(0,Gmax);
    std::cout << "Basis set has " << itsGs.size() << "plane waves G < " << Gmax << std::endl;

    for (List<RVec3>::const_iterator b(itsGs.begin()); b!=itsGs.end(); b++)
        BasisSetImplementation::Insert(new PlaneWaveBF(*b,LatticeVolume));

//
//  Make the integral engine.  Can't do this until all the basis functions and
//  blocks are in place.
//
    if (theMesh)
        TBasisSetImplementation<std::complex<double> >::Insert(new NumericalIE<std::complex<double> >(theMesh));
    else
        TBasisSetImplementation<std::complex<double> >::Insert(new PlaneWaveIE);
};

std::ostream& PlaneWaveBS::Write(std::ostream& os) const
{
    if (!Pretty())
    {
        WriteBasisFunctions(os); //These must be written first.
        BasisSetImplementation::Write(os);
        TBasisSetImplementation<std::complex<double> >::Write(os);
    }
    return os;
}

std::istream& PlaneWaveBS::Read (std::istream& is)
{
    ReadBasisFunctions(is);  //These must be read first.
    BasisSetImplementation::Read(is);
    TBasisSetImplementation<std::complex<double> >::Read(is);
    return is;
}

BasisSet* PlaneWaveBS::Clone() const
{
    return new PlaneWaveBS(*this);
}

BasisSet* PlaneWaveBS::Clone(const RVec3&) const
{
    std::cerr << "Why are you recentering a plane wave basis set!?" << std::endl;
    return new PlaneWaveBS(*this);
}

