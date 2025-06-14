
// stl_io.C Develop stl io code.

#include "gtest/gtest.h"
#include "oml/vector.h"
#include "oml/imp/stream.h"
#include "Base/stl_io.h"
#include "Imp/Cluster/Atom.H"
#include "Imp/Cluster/Molecule.H"
#include <vector>
#include <ostream>
#include <sstream>
#include <fstream>

using std::ostream;
using std::istream;
using std::cout;
using std::endl;
using std::ends;

class STLTesting : public ::testing::TestWithParam<StreamableObject::Mode>
{
public:
    STLTesting() :v(10,3),a(10)
    {
        Fill(a,3);
        for (int i=0;i<6;i++) s.insert(pow(0.5,i));
    }

    std::vector<int> v;
    std::set<double> s;
    Vector<int> a;
};

TEST_F(STLTesting,AsciiIO)
{
    StreamableObject::Mode mode=StreamableObject::SetToAscii();
    std::cout << a << endl;
    std::cout << v << endl;
    std::cout << s << endl;

    {
        std::ostringstream os;
        os << a;
        EXPECT_EQ(os.str(),"1 6VectorIiE 1 10 3 3 3 3 3 3 3 3 3 3 ");
    }
    {
        std::ostringstream os;
        os << v;
        EXPECT_EQ(os.str(),"1 St6vectorIiSaIiEE 10 3 3 3 3 3 3 3 3 3 3 ");
    }
    {
        std::ostringstream os;
        os << s;
        EXPECT_EQ(os.str(),"1 St3setIdSt4lessIdESaIdEE 6 0.03125 0.0625 0.125 0.25 0.5 1 ");
    } 
    StreamableObject::SetOutputMode(mode);    
}


TEST_P(STLTesting,FileIO)
{   
    StreamableObject::SetOutputMode(GetParam());
    //
    //  Stream data to a file
    //
    std::ofstream of("t.dat");
    of << a << v << s;
    of.close();
    //
    //  read it all back in
    //
    Vector<int> a2;
    std::vector<int> v2;
    std::set<double> s2;
    std::ifstream in("t.dat");
    in >> a2 >> v2 >> s2;
    //
    //  Now check in ascii mode.
    //
    StreamableObject::Mode mode=StreamableObject::SetToAscii();
    {
        std::ostringstream os;
        os << a2;
        EXPECT_EQ(os.str(),"1 6VectorIiE 1 10 3 3 3 3 3 3 3 3 3 3 ");
    }
    {
        std::ostringstream os;
        os << v2;
        EXPECT_EQ(os.str(),"1 St6vectorIiSaIiEE 10 3 3 3 3 3 3 3 3 3 3 ");
    }
    { 
        std::ostringstream os;
        os << s2;
        EXPECT_EQ(os.str(),"1 St3setIdSt4lessIdESaIdEE 6 0.03125 0.0625 0.125 0.25 0.5 1 ");
    }
    StreamableObject::SetOutputMode(mode);
}

INSTANTIATE_TEST_CASE_P(FileIO,STLTesting,
                        ::testing::Values(StreamableObject::ascii,StreamableObject::binary));

                        


