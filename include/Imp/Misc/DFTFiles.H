// File: DFTFiles.H  Handle files for the DFT system.
#ifndef _DFTFiles_H_
#define _DFTFiles_H_

#include <string>

class IrrepBasisSet;
class Cluster;
class Mesh;
class Hamiltonian;
class WaveFunction;

//-------------------------------------------------------------------
//
//  Encapsulate management of files names and extensions for the
//  DFT system.
//
class DFTFiles
{
public:
    DFTFiles(const std::string& stub) : itsStub(stub) {};
//
//  Do the files exist.
//
    bool     OBasisSetFile() const;
    bool     CBasisSetFile() const;
    bool     XBasisSetFile() const;
    bool       ClusterFile() const;
    bool          MeshFile() const;
    bool         AtomsFile() const;
    bool ChargeDensityFile() const;
    bool HamiltonianFile() const;
    bool WaveFunctionFile() const;
//
//  Add file name extensions to the stub.
//
    std::string GetDataBaseFN     () const
    {
        return itsStub+".db";
    }
    std::string GetBasisDataFN    () const
    {
        return itsStub+".bsd";
    }
    std::string GetOBasisSetFN    () const
    {
        return itsStub+".obs";
    }
    std::string GetCBasisSetFN    () const
    {
        return itsStub+".cbs";
    }
    std::string GetXBasisSetFN    () const
    {
        return itsStub+".xbs";
    }
    std::string GetClusterFN      () const
    {
        return itsStub+".cl";
    }
    std::string GetAtomsFN        () const
    {
        return itsStub+".atm";
    }
    std::string GetMeshFN         () const
    {
        return itsStub+".msh";
    }
    std::string GetChargeDensityFN() const
    {
        return itsStub+".cd";
    }
    std::string GetHamiltonianFN  () const
    {
        return itsStub+".ham";
    }
    std::string GetWaveFunctionFN () const
    {
        return itsStub+".wav";
    }

private:
    std::string itsStub; //File name stub associated with a given calculation.
};

#endif // _DFTFiles_H_
