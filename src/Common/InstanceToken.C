// File: Common/InstanceToken.C  An object's INSTANCE identity: a never-reused serial, and a weak handle that says whether it still lives.
//
// doc/IntegralCacheOwnershipPlan.md §3b.  A cache that keys a table by ANOTHER object (an orbital block's ⟨ab|c⟩ against a fit basis c)
// needs that object's identity.  Not its pointer: an address is reused after the object is freed, so a new partner at the old address
// would get the dead one's entries (the June 2026 density-pointer bug, M_DFT_Water diverging at iteration 5).  Not a content string:
// the key would be huge (a δ fit's point set) and two equal partners built apart are still two partners.  So a serial from a
// process-wide counter, stamped at construction and never reused, held in a shared_ptr the owner keeps: the cache holds a weak_ptr
// beside the serial and can tell on its own (expired()) that the partner is gone and its entries can never be hit again.
//
// A 64-bit serial, not a GUID: no cache key outlives its process (the stores are RAM only), so uniqueness inside one process is all a
// key needs, and a serial is cheap, ordered and readable in a report.
module;
#include <atomic>
#include <cstdint>
#include <memory>
export module qchem.InstanceToken;

export namespace qchem
{

class InstanceToken
{
public:
    InstanceToken() : itsLife(std::make_shared<const uint64_t>(Next())) {}
    //! A copy is a NEW instance: it mints its own serial (an immutable exact copy sharing the serial would buy nothing, §3b).
    InstanceToken(const InstanceToken&) : InstanceToken() {}
    //! Assignment keeps this object's own identity: the identity belongs to the object, not to its value.
    InstanceToken& operator=(const InstanceToken&) {return *this;}

    uint64_t Serial() const {return *itsLife;}
    //! The weak handle a cache keeps beside the serial: expired exactly when the owning object is gone.
    std::weak_ptr<const uint64_t> Alive() const {return itsLife;}

private:
    static uint64_t Next()
    {
        static std::atomic<uint64_t> theSerial{0};
        return ++theSerial;
    }
    std::shared_ptr<const uint64_t> itsLife;   //!< holds the serial; its control block is what the weak handles watch
};

} // namespace qchem
