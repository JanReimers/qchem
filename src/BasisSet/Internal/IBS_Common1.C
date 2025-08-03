// File: IBS_Common1.C  Irrep Basis set common implementation.
module;
#include <cassert>
#include <iosfwd>
class DiracIntegralTests;

export module qchem.BasisSet.Internal.IBS_Common1;
export import qchem.Irrep_BS;
import qchem.LASolver;
import qchem.BasisSet.Internal.HeapDB;
import qchem.BasisSet.Internal.IEClient;
import qchem.BasisSet.Internal.HeapDB;
import qchem.Fit_IBS;
import qchem.DFT_IBS;
import qchem.HF_IBS;
import qchem.DHF_IBS;

import qchem.BasisSet.Internal.Integrals;

import qchem.Symmetry;
import Common.UniqueIDImp;


