// File: DFTFiles.C  Handle files for the DFT system.

#include "Imp/Misc/DFTFiles.H"
#include "Imp/Misc/Unpickle.H"
#include "Imp/Cluster/UnitCell.H"

#include <Irrep_BS.H>
#include <Hamiltonian.H>
#include <ChargeDensity.H>
#include <Cluster.H>
#include <Mesh.H>
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

