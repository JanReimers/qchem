// qchem.HDF5: the saved-SCF-state carrier (CK-1).  Round trips through a real file, the h5py-facing layout
// facts (C order, the {r,i} complex compound), and the one fallible call (Open on a missing file).
#include <gtest/gtest.h>
#include <complex>
#include <cstdio>
#include <cstdint>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <vector>
import qchem.HDF5;

using namespace qchem;
namespace {

std::string TempPath(const std::string& leaf) {return (std::filesystem::path(testing::TempDir())/leaf).string();}

TEST(HDF5, DatasetsAttributesAndGroupsRoundTrip)
{
    const std::string path=TempPath("qchem_hdf5_roundtrip.h5");
    const std::vector<double>               r{1.0, -2.5, 3.25, 4.0, 5.5, -6.0};
    const std::vector<std::complex<double>> c{{1,2}, {-3,0.5}, {0,-1}, {7,8}};
    const std::vector<std::int64_t>         i{3, -1, 4, 1, 5};
    {
        H5::File f=H5::File::Create(path);
        f.SetAttr("format", "qchem-test 1");
        f.SetAttr("energy", -7.25);
        f.SetAttr("count",  std::int64_t(42));
        f.SetAttr("energy", -8.5);                             // an attribute is a setting: last write wins
        H5::Group b=f.CreateGroup("blocks");
        H5::Group b0=b.CreateGroup("0");
        b0.Write("D", r, {2,3});                               // C order: row 0 = {1,-2.5,3.25}
        b0.Write("C", c, {2,2});
        b0.Write("idx", i);
        f.Flush();
    }
    auto o=H5::File::Open(path);
    ASSERT_TRUE(o) << o.Error();
    H5::File f=o.TakeValue();
    EXPECT_EQ(f.AttrString("format"), "qchem-test 1");
    EXPECT_DOUBLE_EQ(f.AttrReal("energy"), -8.5);
    EXPECT_EQ(f.AttrInt("count"), 42);
    EXPECT_TRUE (f.HasAttr("count"));
    EXPECT_FALSE(f.HasAttr("nope"));
    ASSERT_TRUE(f.Has("blocks"));
    EXPECT_EQ(f.Children(), (std::vector<std::string>{"blocks"}));
    EXPECT_EQ(f.AttrNames(), (std::vector<std::string>{"count","energy","format"}));   // NAME order
    EXPECT_TRUE (f.AttrIsString("format"));
    EXPECT_FALSE(f.AttrIsString("energy"));
    H5::Group b0=f.OpenGroup("blocks").OpenGroup("0");
    EXPECT_EQ(b0.Shape("D"), (std::vector<size_t>{2,3}));
    EXPECT_EQ(b0.Children(), (std::vector<std::string>{"C","D","idx"}));
    EXPECT_FALSE(b0.IsComplex("D"));
    EXPECT_TRUE (b0.IsComplex("C"));
    EXPECT_EQ(b0.ReadReal("D"), r);
    EXPECT_EQ(b0.ReadComplex("C"), c);
    EXPECT_EQ(b0.ReadInt("idx"), i);
    // A real dataset promotes exactly; a complex one refuses to be read as real (never a silent Re()).
    const auto promoted=b0.ReadComplex("D");
    for (size_t k=0;k<r.size();k++) EXPECT_EQ(promoted[k], std::complex<double>(r[k], 0.0));
    EXPECT_THROW(b0.ReadReal("C"), std::runtime_error);
    EXPECT_THROW(b0.ReadReal("missing"), std::runtime_error);
    EXPECT_THROW(f.OpenGroup("missing"), std::runtime_error);
    std::filesystem::remove(path);
}

TEST(HDF5, OpenFailsAsAValueOnAMissingOrForeignFile)
{
    auto missing=H5::File::Open(TempPath("qchem_hdf5_does_not_exist.h5"));
    EXPECT_FALSE(missing);
    const std::string txt=TempPath("qchem_hdf5_not_hdf5.h5");
    { FILE* fp=std::fopen(txt.c_str(), "w"); std::fputs("not an HDF5 file\n", fp); std::fclose(fp); }
    EXPECT_FALSE(H5::File::Open(txt));
    std::filesystem::remove(txt);
}

TEST(HDF5, ShapeMismatchAndDuplicateNamesThrow)
{
    const std::string path=TempPath("qchem_hdf5_errors.h5");
    H5::File f=H5::File::Create(path);
    EXPECT_THROW(f.Write("x", std::vector<double>{1,2,3}, {2,2}), std::runtime_error);
    f.Write("y", std::vector<double>{1,2});
    EXPECT_THROW(f.Write("y", std::vector<double>{1,2}), std::runtime_error);
    f.CreateGroup("g");
    EXPECT_THROW(f.CreateGroup("g"), std::runtime_error);
    std::filesystem::remove(path);
}

} //namespace
