// File: DFTFiles.C  Handle files for the DFT system.

#include "Misc/DFTFiles.H"
#include "Misc/Unpickle.H"
#include "Imp/Containers/ptr_vector.h"
#include "Imp/Containers/ptr_vector_io.h"
#include "Imp/Cluster/UnitCell.H"

#include <BasisSet.H>
#include <Hamiltonian.H>
#include <ChargeDensity.H>
#include <Cluster.H>
#include "Mesh/Mesh.H"
#include <WaveFunction.H>

#include <fstream>
#include <stdlib.h>
#include <iostream>

bool DFTFiles::OBasisSetFile() const
{
    return std::ifstream(GetOBasisSetFN().c_str()).good();
}

bool DFTFiles::CBasisSetFile() const
{
    return std::ifstream(GetCBasisSetFN().c_str()).good();
}

bool DFTFiles::XBasisSetFile() const
{
    return std::ifstream(GetXBasisSetFN().c_str()).good();
}

bool DFTFiles::ClusterFile() const
{
    return std::ifstream(GetClusterFN().c_str()).good();
}

bool DFTFiles::MeshFile() const
{
    return std::ifstream(GetMeshFN().c_str()).good();
}

bool DFTFiles::AtomsFile() const
{
    return std::ifstream(GetAtomsFN().c_str()).good();
}

bool DFTFiles::ChargeDensityFile() const
{
    return std::ifstream(GetChargeDensityFN().c_str()).good();
}


bool DFTFiles::HamiltonianFile() const
{
    return std::ifstream(GetHamiltonianFN().c_str()).good();
}

bool DFTFiles::WaveFunctionFile() const
{
    return std::ifstream(GetWaveFunctionFN().c_str()).good();
}


void DFTFiles::UnpickleOBasisSet(optr_vector1<IrrepBasisSet*>& bs) const
{
    std::cout << "Reading basis set from file " << GetOBasisSetFN() << std::endl;
    std::ifstream in(GetOBasisSetFN().c_str());
    if (!in)
    {
        std::cerr << "Can't open orbital basis set file " << GetOBasisSetFN() << std::endl;
        exit (-1);
    }
    in >> bs;
}

IrrepBasisSet* DFTFiles::UnpickleCBasisSet() const
{
    IrrepBasisSet* ret=0;
    UnPickle(ret,GetCBasisSetFN().c_str(),"Charge density fitting basis set");
    return ret;
}

IrrepBasisSet* DFTFiles::UnpickleXBasisSet() const
{
    IrrepBasisSet* ret=0;
    UnPickle(ret,GetXBasisSetFN().c_str(),"Exchange fitting basis set");
    return ret;
}

Cluster*  DFTFiles::UnpickleCluster() const
{
    Cluster* ret=0;
    UnPickle(ret,GetClusterFN().c_str(),"Cluster");
    return ret;
}

Mesh* DFTFiles::UnpickleMesh() const
{
    Mesh* ret=0;
    UnPickle(ret,GetMeshFN().c_str(),"numerical integration mesh");
    return ret;
}

ChargeDensity* DFTFiles::UnpickleChargeDensity() const
{
    ChargeDensity* ret=0;
    UnPickle(ret,GetChargeDensityFN().c_str(),"Charge density");
    return ret;
}

Hamiltonian* DFTFiles::UnpickleHamiltonian() const
{
    Hamiltonian* ret=0;
    UnPickle(ret,GetHamiltonianFN().c_str(),"Hamiltonian");
    return ret;
}

WaveFunction* DFTFiles::UnpickleWaveFunction() const
{
    WaveFunction* ret=0;
    UnPickle(ret,GetWaveFunctionFN().c_str(),"WaveFunction");
    return ret;
}

