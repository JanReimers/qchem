// File: Common/HDF5.C  A thin RAII carrier over the serial HDF5 C library -- the on-disk format of a saved
// SCF state (doc/OpenWork.md §2 "SCF checkpoint/restart", CK-1) and, by the Viz/GUI plan's choice, of the
// run-report sidecar, so `h5py` reads what we write.
//
// WHAT IT IS NOT.  It knows no Blaze type and no physics: datasets are FLAT buffers plus a shape, in C (row-
// major) order -- h5py's order -- so qcCommon stays the pure-utility leaf it is.  A caller holding a column-
// major Blaze matrix flattens it row by row itself; the layout decision belongs to the one who knows the
// data (the SolidCalculation state layer, src/Calculation/SolidState.C).
//
// COMPLEX is the {r,i} compound -- the layout h5py maps to complex128 natively, so a stored D(k) reads in
// Python as a complex array with no conversion.
//
// ERRORS.  Every failure THROWS (the project's `throw is a marker` policy) with the path and the name, and
// HDF5's own stderr error stack is silenced: the thrown message is the report.  The one exception is
// File::Open, which returns an Outcome -- a missing or foreign file is the ordinary "no saved state yet"
// case a restarting campaign must act on, not a broken invariant.
module;
#include <complex>
#include <cstdint>
#include <string>
#include <vector>
export module qchem.HDF5;
export import qchem.Outcome;

export namespace qchem::H5
{

//! \brief A group (a directory of datasets, sub-groups and attributes).  Move-only; closes itself.
class Group
{
public:
    Group(Group&&) noexcept;
    Group& operator=(Group&&) noexcept;
    Group(const Group&)            = delete;
    Group& operator=(const Group&) = delete;
    virtual ~Group();

    //! \name Sub-groups
    //!@{
    Group CreateGroup(const std::string& name);      //!< THROWS if \a name already exists
    Group OpenGroup  (const std::string& name) const;//!< THROWS if \a name is absent
    bool  Has        (const std::string& name) const;//!< a dataset or a group of that name exists here
    //!@}

    //! \name Datasets: \a data in C order, \f$\prod\f$ \a shape elements (an empty shape = a 1-D array)
    //!@{
    void Write(const std::string& name, const std::vector<double>&               data, std::vector<size_t> shape={});
    void Write(const std::string& name, const std::vector<std::complex<double>>& data, std::vector<size_t> shape={});
    void Write(const std::string& name, const std::vector<std::int64_t>&         data, std::vector<size_t> shape={});
    //! Is the stored dataset the {r,i} complex compound?
    bool IsComplex(const std::string& name) const;
    //! The stored shape, for a caller that must check it before trusting the data.
    std::vector<size_t> Shape(const std::string& name) const;
    std::vector<double>               ReadReal   (const std::string& name) const;  //!< THROWS on a complex dataset
    std::vector<std::complex<double>> ReadComplex(const std::string& name) const;  //!< a REAL dataset is promoted
    std::vector<std::int64_t>         ReadInt    (const std::string& name) const;
    //!@}

    //! \name Attributes: scalars and strings on this group
    //!@{
    void SetAttr(const std::string& name, double);
    void SetAttr(const std::string& name, std::int64_t);
    void SetAttr(const std::string& name, const std::string&);
    void SetAttr(const std::string& name, const char* s) {SetAttr(name, std::string(s));}
    bool         HasAttr   (const std::string& name) const;
    double       AttrReal  (const std::string& name) const;
    std::int64_t AttrInt   (const std::string& name) const;
    std::string  AttrString(const std::string& name) const;
    //!@}

protected:
    Group(std::int64_t id, std::string path);
    std::int64_t itsId = -1;   //!< the hid_t (int64 since HDF5 1.10; the Imp static_asserts it)
    std::string  itsPath;      //!< "file.h5:/blocks/3" -- for the error messages only
};

//! \brief A file, which IS its root group.
class File : public Group
{
public:
    //! Create (TRUNCATING an existing file).  THROWS when the path cannot be written.
    static File Create(const std::string& path);
    //! Open read-only.  FAILS (a value) when the file is absent or is not HDF5 -- the ordinary "nothing saved
    //! yet" case of a restarting campaign.
    static Outcome<File,std::string> Open(const std::string& path);
    File(File&&) noexcept;
    File& operator=(File&&) noexcept;
    ~File() override;
    //! Push buffered writes to disk (a checkpoint writer calls it before renaming the file into place).
    void Flush();
private:
    File(std::int64_t file, std::int64_t root, std::string path);
    std::int64_t itsFile = -1;
};

} //export namespace qchem::H5
