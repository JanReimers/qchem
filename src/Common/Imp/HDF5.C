// File: Common/Imp/HDF5.C  The serial HDF5 C library behind qchem.HDF5.
module;
#include <hdf5.h>
#include <complex>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>
module qchem.HDF5;

static_assert(std::is_same_v<hid_t, std::int64_t>, "qchem.HDF5 stores hid_t as int64 (HDF5 >= 1.10)");

namespace qchem::H5
{

//---------------------------------------------------------------------------------------------------
//  Plumbing
//---------------------------------------------------------------------------------------------------
namespace
{
// HDF5 prints its whole error stack to stderr on every failed call -- including the EXPECTED ones (probing a
// file that is not there).  Silence it once; every failure below throws with its own message instead.
struct Silence { Silence() {H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr);} };
void Quiet() {static Silence s;}

[[noreturn]] void Fail(const std::string& path, const std::string& what)
{
    throw std::runtime_error("HDF5 ("+path+"): "+what);
}
hid_t Check(hid_t id, const std::string& path, const std::string& what)
{
    if (id<0) Fail(path, what);
    return id;
}
void Check(herr_t e, const std::string& path, const std::string& what, int)
{
    if (e<0) Fail(path, what);
}

//! One owned handle of any kind, closed by the matching H5?close.
class Handle
{
public:
    using closer_t = herr_t(*)(hid_t);
    Handle(hid_t id, closer_t c) : itsId(id), itsClose(c) {}
    ~Handle() {if (itsId>=0) itsClose(itsId);}
    Handle(const Handle&) = delete;
    Handle& operator=(const Handle&) = delete;
    operator hid_t() const {return itsId;}
private:
    hid_t    itsId;
    closer_t itsClose;
};

//! The {r,i} compound: h5py reads it as complex128.
hid_t ComplexType()
{
    hid_t t=H5Tcreate(H5T_COMPOUND, sizeof(std::complex<double>));
    H5Tinsert(t, "r", 0,              H5T_NATIVE_DOUBLE);
    H5Tinsert(t, "i", sizeof(double), H5T_NATIVE_DOUBLE);
    return t;
}
hid_t StringType()
{
    hid_t t=H5Tcopy(H5T_C_S1);
    H5Tset_size(t, H5T_VARIABLE);
    H5Tset_cset(t, H5T_CSET_UTF8);
    return t;
}

size_t Count(const std::vector<size_t>& shape)
{
    size_t n=1;
    for (size_t d : shape) n*=d;
    return n;
}

template <class T> void WriteDataset(hid_t g, const std::string& path, const std::string& name,
                                     const std::vector<T>& data, std::vector<size_t> shape, hid_t memType)
{
    if (shape.empty()) shape={data.size()};
    if (Count(shape)!=data.size())
        Fail(path, "dataset '"+name+"': "+std::to_string(data.size())+" values do not fill the stated shape");
    std::vector<hsize_t> dims(shape.begin(), shape.end());
    Handle space(Check(H5Screate_simple(int(dims.size()), dims.data(), nullptr), path, "dataspace for '"+name+"'"), H5Sclose);
    Handle ds(Check(H5Dcreate2(g, name.c_str(), memType, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT),
                    path, "cannot create dataset '"+name+"' (does it already exist?)"), H5Dclose);
    if (!data.empty())
        Check(H5Dwrite(ds, memType, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data()), path, "writing '"+name+"'", 0);
}

std::vector<size_t> DatasetShape(hid_t ds, const std::string& path, const std::string& name)
{
    Handle space(Check(H5Dget_space(ds), path, "dataspace of '"+name+"'"), H5Sclose);
    const int rank=H5Sget_simple_extent_ndims(space);
    if (rank<0) Fail(path, "rank of '"+name+"'");
    std::vector<hsize_t> dims(rank);
    if (rank>0) H5Sget_simple_extent_dims(space, dims.data(), nullptr);
    return std::vector<size_t>(dims.begin(), dims.end());
}

template <class T> std::vector<T> ReadDataset(hid_t g, const std::string& path, const std::string& name, hid_t memType)
{
    Handle ds(Check(H5Dopen2(g, name.c_str(), H5P_DEFAULT), path, "no dataset '"+name+"'"), H5Dclose);
    std::vector<T> out(Count(DatasetShape(ds, path, name)));
    if (!out.empty())
        Check(H5Dread(ds, memType, H5S_ALL, H5S_ALL, H5P_DEFAULT, out.data()), path, "reading '"+name+"'", 0);
    return out;
}

void WriteAttr(hid_t g, const std::string& path, const std::string& name, hid_t type, const void* value)
{
    if (H5Aexists(g, name.c_str())>0) H5Adelete(g, name.c_str());   // an attribute is a SETTING: last write wins
    Handle space(Check(H5Screate(H5S_SCALAR), path, "scalar dataspace"), H5Sclose);
    Handle a(Check(H5Acreate2(g, name.c_str(), type, space, H5P_DEFAULT, H5P_DEFAULT), path,
                   "cannot create attribute '"+name+"'"), H5Aclose);
    Check(H5Awrite(a, type, value), path, "writing attribute '"+name+"'", 0);
}
void ReadAttr(hid_t g, const std::string& path, const std::string& name, hid_t type, void* value)
{
    Handle a(Check(H5Aopen(g, name.c_str(), H5P_DEFAULT), path, "no attribute '"+name+"'"), H5Aclose);
    Check(H5Aread(a, type, value), path, "reading attribute '"+name+"'", 0);
}
} //anonymous namespace

//---------------------------------------------------------------------------------------------------
//  Group
//---------------------------------------------------------------------------------------------------
Group::Group(std::int64_t id, std::string path) : itsId(id), itsPath(std::move(path)) {}
Group::Group(Group&& o) noexcept : itsId(std::exchange(o.itsId, -1)), itsPath(std::move(o.itsPath)) {}
Group& Group::operator=(Group&& o) noexcept
{
    if (this!=&o)
    {
        if (itsId>=0) H5Gclose(itsId);
        itsId  =std::exchange(o.itsId, -1);
        itsPath=std::move(o.itsPath);
    }
    return *this;
}
Group::~Group() {if (itsId>=0) H5Gclose(itsId);}

static std::string Child(const std::string& path, const std::string& name)
{
    return path + (path.back()=='/' ? "" : "/") + name;
}
Group Group::CreateGroup(const std::string& name)
{
    hid_t g=Check(H5Gcreate2(itsId, name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT), itsPath,
                  "cannot create group '"+name+"' (does it already exist?)");
    return Group(g, Child(itsPath, name));
}
Group Group::OpenGroup(const std::string& name) const
{
    hid_t g=Check(H5Gopen2(itsId, name.c_str(), H5P_DEFAULT), itsPath, "no group '"+name+"'");
    return Group(g, Child(itsPath, name));
}
bool Group::Has(const std::string& name) const {return H5Lexists(itsId, name.c_str(), H5P_DEFAULT)>0;}

static herr_t CollectLink(hid_t, const char* name, const H5L_info2_t*, void* out)
{
    static_cast<std::vector<std::string>*>(out)->push_back(name);
    return 0;
}
std::vector<std::string> Group::Children() const
{
    std::vector<std::string> out;
    Check(H5Literate2(itsId, H5_INDEX_NAME, H5_ITER_INC, nullptr, CollectLink, &out), itsPath, "listing the group", 0);
    return out;
}

void Group::Write(const std::string& name, const std::vector<double>& data, std::vector<size_t> shape)
{
    WriteDataset(itsId, itsPath, name, data, std::move(shape), H5T_NATIVE_DOUBLE);
}
void Group::Write(const std::string& name, const std::vector<std::complex<double>>& data, std::vector<size_t> shape)
{
    Handle t(ComplexType(), H5Tclose);
    WriteDataset(itsId, itsPath, name, data, std::move(shape), t);
}
void Group::Write(const std::string& name, const std::vector<std::int64_t>& data, std::vector<size_t> shape)
{
    WriteDataset(itsId, itsPath, name, data, std::move(shape), H5T_NATIVE_INT64);
}

bool Group::IsComplex(const std::string& name) const
{
    Handle ds(Check(H5Dopen2(itsId, name.c_str(), H5P_DEFAULT), itsPath, "no dataset '"+name+"'"), H5Dclose);
    Handle t(Check(H5Dget_type(ds), itsPath, "type of '"+name+"'"), H5Tclose);
    return H5Tget_class(t)==H5T_COMPOUND;
}
std::vector<size_t> Group::Shape(const std::string& name) const
{
    Handle ds(Check(H5Dopen2(itsId, name.c_str(), H5P_DEFAULT), itsPath, "no dataset '"+name+"'"), H5Dclose);
    return DatasetShape(ds, itsPath, name);
}
std::vector<double> Group::ReadReal(const std::string& name) const
{
    if (IsComplex(name)) Fail(itsPath, "dataset '"+name+"' is COMPLEX; read it with ReadComplex");
    return ReadDataset<double>(itsId, itsPath, name, H5T_NATIVE_DOUBLE);
}
std::vector<std::complex<double>> Group::ReadComplex(const std::string& name) const
{
    if (!IsComplex(name))                                       // a real dataset promotes exactly
    {
        const std::vector<double> r=ReadReal(name);
        return std::vector<std::complex<double>>(r.begin(), r.end());
    }
    Handle t(ComplexType(), H5Tclose);
    return ReadDataset<std::complex<double>>(itsId, itsPath, name, t);
}
std::vector<std::int64_t> Group::ReadInt(const std::string& name) const
{
    return ReadDataset<std::int64_t>(itsId, itsPath, name, H5T_NATIVE_INT64);
}

void Group::SetAttr(const std::string& name, double v)       {WriteAttr(itsId, itsPath, name, H5T_NATIVE_DOUBLE, &v);}
void Group::SetAttr(const std::string& name, std::int64_t v) {WriteAttr(itsId, itsPath, name, H5T_NATIVE_INT64,  &v);}
void Group::SetAttr(const std::string& name, const std::string& v)
{
    Handle t(StringType(), H5Tclose);
    const char* s=v.c_str();
    WriteAttr(itsId, itsPath, name, t, &s);
}
bool Group::HasAttr(const std::string& name) const {return H5Aexists(itsId, name.c_str())>0;}
bool Group::AttrIsString(const std::string& name) const
{
    Handle a(Check(H5Aopen(itsId, name.c_str(), H5P_DEFAULT), itsPath, "no attribute '"+name+"'"), H5Aclose);
    Handle t(Check(H5Aget_type(a), itsPath, "type of attribute '"+name+"'"), H5Tclose);
    return H5Tget_class(t)==H5T_STRING;
}
static herr_t CollectAttr(hid_t, const char* name, const H5A_info_t*, void* out)
{
    static_cast<std::vector<std::string>*>(out)->push_back(name);
    return 0;
}
std::vector<std::string> Group::AttrNames() const
{
    std::vector<std::string> out;
    Check(H5Aiterate2(itsId, H5_INDEX_NAME, H5_ITER_INC, nullptr, CollectAttr, &out), itsPath, "listing the attributes", 0);
    return out;
}
double Group::AttrReal(const std::string& name) const
{
    double v=0.0;
    ReadAttr(itsId, itsPath, name, H5T_NATIVE_DOUBLE, &v);
    return v;
}
std::int64_t Group::AttrInt(const std::string& name) const
{
    std::int64_t v=0;
    ReadAttr(itsId, itsPath, name, H5T_NATIVE_INT64, &v);
    return v;
}
std::string Group::AttrString(const std::string& name) const
{
    Handle t(StringType(), H5Tclose);
    char* s=nullptr;
    ReadAttr(itsId, itsPath, name, t, &s);
    std::string out = s ? std::string(s) : std::string();
    if (s) H5free_memory(s);
    return out;
}

//---------------------------------------------------------------------------------------------------
//  File
//---------------------------------------------------------------------------------------------------
File::File(std::int64_t file, std::int64_t root, std::string path) : Group(root, path+":/"), itsFile(file) {}
File::File(File&& o) noexcept : Group(std::move(o)), itsFile(std::exchange(o.itsFile, -1)) {}
File& File::operator=(File&& o) noexcept
{
    if (this!=&o)
    {
        Group::operator=(std::move(o));
        if (itsFile>=0) H5Fclose(itsFile);
        itsFile=std::exchange(o.itsFile, -1);
    }
    return *this;
}
File::~File()
{
    if (itsId>=0) {H5Gclose(itsId); itsId=-1;}   // the root group first, so the file closes for real
    if (itsFile>=0) H5Fclose(itsFile);
}

File File::Create(const std::string& path)
{
    Quiet();
    hid_t f=Check(H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT), path, "cannot create the file");
    hid_t g=H5Gopen2(f, "/", H5P_DEFAULT);
    if (g<0) {H5Fclose(f); Fail(path, "cannot open the root group");}
    return File(f, g, path);
}

Outcome<File,std::string> File::Open(const std::string& path)
{
    using O=Outcome<File,std::string>;
    Quiet();
    if (H5Fis_accessible(path.c_str(), H5P_DEFAULT)<=0)
        return O::Fail("'"+path+"' is absent or is not an HDF5 file");
    hid_t f=H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
    if (f<0) return O::Fail("'"+path+"' could not be opened read-only");
    hid_t g=H5Gopen2(f, "/", H5P_DEFAULT);
    if (g<0) {H5Fclose(f); return O::Fail("'"+path+"' has no root group");}
    return O::Ok(File(f, g, path));
}

void File::Flush() {Check(H5Fflush(itsFile, H5F_SCOPE_GLOBAL), itsPath, "flush", 0);}

} //namespace qchem::H5
