// File: UnitTests/A_Cache4.C  Unit test the Atom Cach4 system for storing 4 index redial Slater integrals.
#include "gtest/gtest.h"
#include <gmock/gmock.h>
#include <iostream>
#include <map>
#include <string>
#include <memory>
#include "../src/forward.H"
// #include <blaze/Math.h>
// #include <nlohmann/json.hpp>
using std::cout;
using std::endl;

import qchem.BasisSet1.DB_Cache;
import qchem.BasisSet1.Atom.Evaluators.IBS; 
import qchem.BasisSet1.Atom.Evaluators.Gaussian.IBS; 
import qchem.BasisSet1.Atom.Evaluators.Internal.Rk;
import qchem.BasisSet1.Atom.Evaluators.Internal.ExponentGrouper;
import qchem.BasisSet1.Atom.Evaluators.Gaussian.Internal.Rk; 
import qchem.stl_io;
// import qchem.BasisSet1.Internal.Cache4;



//
//  Cache object based on four unsigned integer indices.  Client code using the caching
//  should derive from this class and overload all the protected functions with function forwarding.
//  Use covariant return types for the loop_4 overload. 
//  Derived class also needs to supply a Create function.
//
class Cache41
{
public: 
    virtual ~Cache41() {};
    virtual void Register(IBS_Evaluator * eval)=0;

    // void       loop_1(size_t i1) const;
    // void       loop_2(size_t i2) const;
    // void       loop_3(size_t i3) const;
    // virtual const Cacheable* loop_4(size_t i4) const;
    
    
private:
    virtual const Cacheable* Create(size_t i1,size_t i2,size_t i3,size_t i4) const=0;

    typedef std::map<size_t,const Cacheable*> cache_4; 
    typedef std::map<size_t,cache_4> cache_3; 
    typedef std::map<size_t,cache_3> cache_2; 
    typedef std::map<size_t,cache_2> cache_t; 
    
    mutable cache_t cache;
    mutable cache_2* i1_cache;
    mutable cache_3* i2_cache;
    mutable cache_4* i3_cache;
    mutable size_t i1,i2,i3,i4; //Current indexes
};

class Gaussian_Cache4 : public  Cache41
{
public:
    using IBS_Evaluator_t = Gaussian_IBS_Evaluator;
    virtual void Register(IBS_Evaluator * eval)
    {
        assert(eval);
        eval->Register(&grouper);
    }
    virtual Rk*  Create (size_t ia,size_t ic,size_t ib,size_t id) const
    {
        return new Gaussian::RkEngine(
            grouper.unique_esv[ia]+grouper.unique_esv[ib],
            grouper.unique_esv[ic]+grouper.unique_esv[id],
            grouper.LMax(ia,ib,ic,id));
    }
private:
    friend class Cache4Tests;
    ExponentGrouper grouper;
};


class Cache4_DB
{
    using key_t=std::string;
    using val_t=std::unique_ptr<Cache41>;
public:
    // TODO: Evaluator shoudl support Type() instead of Name()
    // TODO: use the evaluator to create this eval->CreatCache4();
    Cache41* Factory(const key_t& types)
    {
        Cache41* ret=0;
        if (types=="Spherical Gaussian ")
            ret=new Gaussian_Cache4();

        assert(ret);
        return ret;
    }

    void Register(IBS_Evaluator* eval)
    {
        key_t key=eval->Name(); //This is really a type
        auto it=itsCache4s.find(key);
        if (it==itsCache4s.end())
        {
            const auto [iterator, success]=itsCache4s.insert({key,val_t(Factory(key))});
            assert(success);
            it=iterator;
        }
        it->second->Register(eval);
    }

    const Cache41* Get(const key_t& type) const
    {
        auto it=itsCache4s.find(type);
        assert(it!=itsCache4s.end());
        return it->second.get();
    }

private:
    friend class Cache4Tests;
    std::map<key_t,val_t> itsCache4s;
};

class Cache4Tests : public ::testing::Test
{
public:
    Cache4Tests() 
    {
        auto cache=BasisSet1::theGlobalCache;
        if (cache==0)
            cache=new BasisSet1::IntegralsCache_RAM<double>(true);     

    }

    const ExponentGrouper& GetGrouper(const Gaussian_Cache4* gc)
    {
        assert(gc);
        return gc->grouper; //friend access.
    }

    const std::vector<size_t>& maxls(const Gaussian_Cache4* gc)
    {
        return GetGrouper(gc).maxls;
    }

    const std::vector<size_t>& es_indices(const IBS_Evaluator* e)
    {
        return e->es_indices;
    }

    Cache4_DB itsCache4;
};


TEST_F(Cache4Tests,test1)
{

    auto s=new Gaussian_IBS_Evaluator({1,2,4,8},0);
    auto p=new Gaussian_IBS_Evaluator({.5,1,2.0,4},1);
    auto d=new Gaussian_IBS_Evaluator({.25,.5,1.0,2},2);
    auto f=new Gaussian_IBS_Evaluator({.25,.5,1.0},3);

    itsCache4.Register(s);
    itsCache4.Register(p);
    itsCache4.Register(d);
    itsCache4.Register(f);

    const Cache41* sc=itsCache4.Get(s->Name());
    EXPECT_EQ(sc,itsCache4.Get(p->Name()));
    EXPECT_EQ(sc,itsCache4.Get(d->Name()));
    EXPECT_EQ(sc,itsCache4.Get(f->Name()));

    const Gaussian_Cache4* sgc=dynamic_cast<const Gaussian_Cache4*>(sc);
    EXPECT_TRUE(sgc);

    const ExponentGrouper& gr=GetGrouper(sgc);

    EXPECT_THAT(gr.unique_esv,::testing::ElementsAre(1,2,4,8,.5,.25));
    std::cout << "Exponents: " << gr.unique_esv << std::endl;
    std::cout << "maxls: " << maxls(sgc) << std::endl;
    std::cout << "s->es_indices: " << es_indices(s) << std::endl;
    EXPECT_THAT(gr.unique_esv,::testing::ElementsAre(1,2,4,8,.5,.25));
    EXPECT_THAT(maxls(sgc),::testing::ElementsAre(3,2,1,0,3,3));
    EXPECT_THAT(es_indices(s),::testing::ElementsAre(0,1,2,3));
    EXPECT_THAT(es_indices(p),::testing::ElementsAre(4,0,1,2));
    EXPECT_THAT(es_indices(d),::testing::ElementsAre(5,4,0,1));
    EXPECT_THAT(es_indices(f),::testing::ElementsAre(5,4,0));

}
