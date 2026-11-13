// File: Common/tests/SlotMap.C  The keyed cache's three promises (doc/IntegralCacheOwnershipPlan.md S1): a threaded fill is the
// serial one, a fill that asks another slot does not deadlock, and a density-dependent (Refresh) slot is never filled lazily.
#include "gtest/gtest.h"
#include <atomic>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>
import qchem.SlotMap;

using namespace qchem;

namespace
{
using OnceMap   =SlotMap<int,std::vector<double>,SlotDiscipline::FillOnce>;
using RefreshMap=SlotMap<int,std::vector<double>,SlotDiscipline::Refresh>;
std::vector<double> Build(int k) {std::vector<double> v(100+k); for (size_t i=0;i<v.size();++i) v[i]=k*1000.0+i; return v;}
}

// A Refresh map has no lazy fill AT ALL: the type says so, so the misuse is a compile error, not a run-time check.
template <class M> concept HasLazyFill=requires(const M& m) {m.GetOrFill(0, []{return std::vector<double>{};});};
template <class M> concept HasRefresh =requires(const M& m) {m.BeginRefresh();};
static_assert( HasLazyFill<OnceMap>    && !HasRefresh<OnceMap>);
static_assert(!HasLazyFill<RefreshMap> &&  HasRefresh<RefreshMap>);

TEST(SlotMap, AThreadedFillIsTheSerialOneAndEachSlotIsBuiltOnce)
{
    constexpr int nKeys=16, nThreads=8, nRounds=50;
    OnceMap serial, threaded;
    for (int k=0;k<nKeys;++k) serial.GetOrFill(k, [k]{return Build(k);});

    std::atomic<int> builds{0};
    std::vector<std::thread> pool;
    for (int t=0;t<nThreads;++t)
        pool.emplace_back([&,t]
        {
            for (int r=0;r<nRounds;++r)
            {
                const int k=(t+r)%nKeys;   // every thread asks every key, many times, in a different order
                threaded.GetOrFill(k, [&,k]{builds++; return Build(k);});
            }
        });
    for (auto& th : pool) th.join();

    EXPECT_EQ(builds.load(), nKeys) << "each slot is built ONCE however many threads ask it";
    for (int k=0;k<nKeys;++k)
    {
        ASSERT_NE(threaded.Read(k), nullptr);
        EXPECT_EQ(*threaded.Read(k), *serial.Read(k)) << "key " << k;
    }
    EXPECT_EQ(threaded.Filled(), size_t(nKeys));
    EXPECT_EQ(threaded.Bytes([](const std::vector<double>& v){return v.size()*sizeof(double);}),
              serial  .Bytes([](const std::vector<double>& v){return v.size()*sizeof(double);}));
}

TEST(SlotMap, AFillThatAsksAnotherSlotDoesNotDeadlock)
{
    OnceMap m;
    // slot 2 is built from slot 1, which is built from slot 0: three nested fills of one map
    const auto& v2=m.GetOrFill(2, [&]
    {
        const auto& v1=m.GetOrFill(1, [&]{auto v=m.GetOrFill(0, []{return Build(0);}); v.push_back(1); return v;});
        auto v=v1; v.push_back(2); return v;
    });
    EXPECT_EQ(v2.size(), Build(0).size()+2);
    EXPECT_EQ(m.Filled(), 3u);
}

TEST(SlotMap, AFillThatThrowsLeavesTheSlotUnfilledAndTheNextAskRetries)
{
    OnceMap m;
    EXPECT_THROW(m.GetOrFill(7, []()->std::vector<double>{throw std::runtime_error("no");}), std::runtime_error);
    EXPECT_EQ(m.Read(7), nullptr);
    EXPECT_EQ(m.GetOrFill(7, []{return Build(7);}), Build(7));
}

TEST(SlotMap, ARefreshSlotIsFilledOnlyInsideItsRefreshAndTheReadIsPure)
{
    RefreshMap m;
    EXPECT_EQ(m.Read(1), nullptr) << "a read never builds";
    {
        auto r=m.BeginRefresh();
        EXPECT_THROW((void)m.BeginRefresh(), std::logic_error) << "one refresh at a time";
        int builds=0;
        r.Fill(1, [&]{builds++; return Build(1);});
        r.Fill(1, [&]{builds++; return Build(99);});   // a second key folded onto the same slot: built once
        EXPECT_EQ(builds, 1);
    }
    ASSERT_NE(m.Read(1), nullptr);
    EXPECT_EQ(*m.Read(1), Build(1));
    {
        auto r=m.BeginRefresh();                        // the next density: the old blocks are gone
        EXPECT_EQ(m.Read(1), nullptr);
        r.Fill(2, []{return Build(2);});
    }
    EXPECT_EQ(m.Read(1), nullptr);
    EXPECT_EQ(m.Filled(), 1u);
}

TEST(SlotMap, ACopyHoldsTheFilledValues)
{
    OnceMap m;
    m.GetOrFill(3, []{return Build(3);});
    const OnceMap c(m);
    ASSERT_NE(c.Read(3), nullptr);
    EXPECT_EQ(*c.Read(3), Build(3));
    int builds=0;
    c.GetOrFill(3, [&]{builds++; return Build(3);});
    EXPECT_EQ(builds, 0) << "the copy's slot is already filled";
}
