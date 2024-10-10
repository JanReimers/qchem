// File: Persistance.C  Test the DFT persistance classes

#include "gtest/gtest.h"
#include "Misc/pmstream.h"
#include "Misc/Persistent/PerRef.H"
#include "Misc/UniqueID/UniqueID.H"
#include "Imp/Containers/ptr_vector.h"
#include "Imp/Containers/ptr_vector_io.h"
#include <iostream> 
#include <fstream>

using std::cout;
using std::endl;

//----------------------------------------------------------------------------------------------
//
//  Polymorphic class structure to run tests on
//
class TestBase : public virtual PMStreamableObject
    , public UniqueID
{
public:
    virtual ~TestBase() {};
    virtual std::ostream& Write(std::ostream&) const=0;
    virtual std::istream& Read (std::istream&)      =0;
    virtual bool operator==(const TestBase&) const=0;
    virtual TestBase* Clone() const=0;
    static TestBase* Factory(std::istream&);
};

class D1 : public virtual TestBase
{
public:
    D1(int i) : D1Data(i) {};
    D1() : D1Data(-1) {};
//    ~D1() {cout << "Deleteing D1 " << (void*)this << endl;}
    virtual std::ostream& Write(std::ostream& os) const
    {
        UniqueID::Write(os);
        return os << D1Data << " ";
    }
    virtual std::istream& Read (std::istream& is)
    {
        UniqueID::Read(is);
        is >> D1Data;
//        cout << "D1::Read D1Data=" << D1Data << endl;
        return is;
    }

    virtual bool operator==(const TestBase& other) const
    {
        const D1* cast_other=dynamic_cast<const D1*>(&other);
        return cast_other && D1Data==cast_other->D1Data;
    }

    virtual TestBase* Clone() const {return new D1(*this);}

    int D1Data;
};

class D2 : public D1
{
public:
    D2(int i,const char* text) : D1(i), D2Data(text) {};
    D2() : D2Data("Not_Initialized") {};
//    ~D2() {cout << "Deleteing D2 " << (void*)this << endl;}
    virtual std::ostream& Write(std::ostream& os) const
    {
        os << D2Data << " ";
        return D1::Write(os);
    }
    virtual std::istream& Read (std::istream& is)
    {
        is >> D2Data;
        return D1::Read(is);
    }
    virtual bool operator==(const TestBase& other) const
    {
        const D2* cast_other=dynamic_cast<const D2*>(&other);
        return cast_other && D1::operator==(other) && D2Data==cast_other->D2Data;
    }
    virtual TestBase* Clone() const {return new D2(*this);}
    std::string D2Data;
};


class Owner
    : public virtual TestBase
{
public:
    Owner () : itsRef1() {};
    Owner(TestBase* tb) :  itsRef1(tb) {};
    virtual std::ostream& Write(std::ostream& os) const
    {
        return os << itsRef1;
    }
    virtual std::istream& Read (std::istream& is)
    {
        return is >> itsRef1;
    }
    virtual bool operator==(const TestBase& other) const
    {
        const Owner* cast_other=dynamic_cast<const Owner*>(&other);
        return cast_other && *itsRef1==*cast_other->itsRef1;
    }
    virtual TestBase* Clone() const {return new Owner(*this);}
    PerRef<TestBase> itsRef1;

};

#define TYPE_STRING "TestBase"
#define TYPE TestBase

#include "DFTDataBase/Instance.Ci"

TestBase* TestBase::Factory(std::istream& is)
{
    TestBase* ret=0;
    std::string Name=StreamableObject::PeekAtName(is);
    if (Name==typeid(D1).name()) ret= new D1;
    else if (Name==typeid(D2).name()) ret= new D2;
    else if (Name==typeid(Owner).name()) ret= new Owner;

    if (!ret) std::cerr << "TestBase::Factory: Unknown object type :" << Name << endl;
    return ret;
}



//----------------------------------------------------------------------------------------------
//
//  Testing class
//
class PersistanceTests : public ::testing::Test
{
public:
    PersistanceTests() : file_name("streamtest1.dat") {StreamableObject::SetToBinary();};
    std::string file_name;

    void OutputD1(int i)
    {
        PMStreamableObject* d1=new D1(i);
        std::ofstream out(file_name.c_str());
        out << *d1;
        EXPECT_TRUE(out);
    }

    void OutputD2(int i, const char* text)
    {
        PMStreamableObject* d2=new D2(i,text);
        std::ofstream out(file_name.c_str());
        out << *d2;
        EXPECT_TRUE(out);
    }
};

TEST_F(PersistanceTests, OutputPerRefD1)
{
    PerDB<TestBase> theDB;
    TestBase* d1=new D1(101);
    PerRef<TestBase> pd1(d1);
    std::ofstream out(file_name.c_str());
    out << pd1;
    EXPECT_TRUE(out);
}

TEST_F(PersistanceTests, OutputOwnerD1)
{
    PerDB<TestBase> theDB;
    TestBase* d1=new D1(101);
    Owner o(d1);
    std::ofstream out(file_name.c_str());
    out << o;
    EXPECT_TRUE(out);
}

TEST_F(PersistanceTests, OutputPerDB)
{
    PerDB<TestBase> theDB;
    TestBase* d1=new D1(101);
    Owner o(d1);
    std::ofstream out(file_name.c_str());
    out << theDB;
    EXPECT_TRUE(out);
}

TEST_F(PersistanceTests, OutputAndInputPerDB)
{
    StreamableObject::SetToBinary();
    {
        PerDB<TestBase> theDB;
        TestBase* d1=new D1(101);
        Owner o(d1);
        std::ofstream out(file_name.c_str());
        out << theDB;
        EXPECT_TRUE(out);
        out.close();
    }
    {
        PerDB<TestBase> theDB;
        std::ifstream in(file_name.c_str());
        in >> theDB;
        EXPECT_TRUE(in);
        in.close();
    }
}

TEST_F(PersistanceTests, OutputAndInputPerDBAndOwner)
{
    {
        PerDB<TestBase> theDB;
        TestBase* d1=new D1(101);
        Owner o(d1);
        std::ofstream out(file_name.c_str());
        out << theDB;
        out << o;
        EXPECT_TRUE(out);

    }
    {
        PerDB<TestBase> theDB;
        std::ifstream in(file_name.c_str());
        in >> theDB;
        Owner o;
        in >> o;
    }
}

TEST_F(PersistanceTests, FixUpPointers)
{
    TestBase* d_compare=new D1(101);
    {
        TestBase* d1=d_compare->Clone();
        PerDB<TestBase> theDB;
        Owner o1(d1);

        std::ofstream out(file_name.c_str());
        out << theDB;
        out << o1;
        EXPECT_TRUE(out);
    }
    //theDB.Clear();
    {
        std::ifstream in(file_name.c_str());
        PerDB<TestBase> theDB;
        in >> theDB;
        Owner o;
        in >> o;
        TestBase* d11=o.itsRef1;
        EXPECT_NE(d11,(void*)0);
        EXPECT_TRUE(*d11==*d_compare    );
        D1* d111=dynamic_cast<D1*>(d11);
        EXPECT_NE(d111,(void*)0);
        EXPECT_TRUE(*d111==*d_compare    );

    }
    delete d_compare;

}

TEST_F(PersistanceTests, FixUpPointers2)
{
    StreamableObject::SetToAscii(); //TODO fails in binary mode.

    int N=10 ;
    optr_vector1<TestBase*> v1;
    for (int i=0;i<N;i++)
    {
        v1.push_back(new D1(200+i));
        v1.push_back(new D2(300+i,"MER"));
    }
    {
        PerDB<TestBase> theDB;
        optr_vector1<Owner*> vo;
        for (auto i:v1) vo.push_back(new Owner(i->Clone()));

        std::ofstream out(file_name.c_str());
        out << theDB;
        out << vo;
        EXPECT_TRUE(out);
    }
    //theDB.Clear();
    {
        std::ifstream in(file_name.c_str());
        PerDB<TestBase> theDB;
        in >> theDB;
        optr_vector1<TestBase*> vo1;
        in >> vo1;

        auto iv1=v1.begin();
        for (auto io:vo1)
        {
            const Owner* o=dynamic_cast<const Owner*>(io);
            const TestBase* tb=o->itsRef1;
            EXPECT_NE(tb,(void*)0);
            EXPECT_TRUE(*tb==**iv1);
            iv1++;
        }

    }

}

