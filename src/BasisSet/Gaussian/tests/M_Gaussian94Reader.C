// File: M_Gaussian94Reader.C  The Gaussian-94 basis reader names an element its file does not hold.
//
// Found 2026-10-21: TiO2 on VALENCE_LOWQ_VA (a file with no Ti).  FindAtom ran to EOF and reported "not found" through a bool every
// caller ignored, so the next ReadNext parsed whatever bytes it stood on and died as "'\xef' is not an angular-momentum label".
// FindAtom now THROWS, naming the element and the file.  The file here is written by the test, so no data-file edit can move it.

#include "gtest/gtest.h"
#include <filesystem>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <string>

import qchem.BasisSet.Gaussian.Point.Readers.Gaussian94;
import qchem.BasisSet.Gaussian.Evaluators.PG_Cart_MnD.GaussianRF;
import qchem.Structure;
import qchem.Types;

using namespace qchem;
using namespace qchem::BasisSet::Gaussian;
using namespace qchem::BasisSet::Gaussian::Evaluators::PG_Cart_MnD;

namespace
{
//! A one-element (H) Gaussian-94 file: the "!" title block the reader seeks, one s shell, the "****" terminator.
std::filesystem::path HydrogenOnlyFile()
{
    const auto p=std::filesystem::temp_directory_path() / "qchem_UT_Gaussian94Reader_H_only.bsd";
    std::ofstream f(p);
    f << "header text the reader skips\n!\nBASIS=\"UT: hydrogen only\"\n H   0\n S   1  1.00\n      1.0   1.0\n****\n";
    return p;
}
}

TEST(Gaussian94Reader, AnElementThatIsInTheFileIsFoundAndRead)
{
    Gaussian94Reader reader(HydrogenOnlyFile().string());
    const Atom H(1, 0.0, Vector3D(0,0,0));
    reader.FindAtom(H);
    std::unique_ptr<GaussianRF> rf{reader.ReadNext(H)};
    EXPECT_TRUE(rf) << "the H s shell";
}

TEST(Gaussian94Reader, AnElementThatIsNotInTheFileThrowsNamingTheElementAndTheFile)
{
    const auto file=HydrogenOnlyFile();
    Gaussian94Reader reader(file.string());
    const Atom Ti(22, 0.0, Vector3D(0,0,0));
    try
    {
        reader.FindAtom(Ti);
        FAIL() << "FindAtom must throw for an element the file does not hold";
    }
    catch (const std::runtime_error& e)
    {
        const std::string what=e.what();
        EXPECT_NE(what.find("Ti"), std::string::npos) << what;
        EXPECT_NE(what.find(file.string()), std::string::npos) << what;
    }
}
