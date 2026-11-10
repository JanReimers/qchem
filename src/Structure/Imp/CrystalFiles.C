// File: Structure/Imp/CrystalFiles.C  CIF / POSCAR / extended xyz <-> one materials.json entry (the translators).
module;
#include <algorithm>
#include <array>
#include <cctype>
#include <chrono>
#include <cmath>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <nlohmann/json.hpp>
module qchem.CrystalFiles;
import qchem.PeriodicTable;
import qchem.Matrix3D;
import qchem.Types;
import qchem.UnitCell;
import qchem.StructureData;                      // CellFromEntry: the one parser every entry (an imported one too) goes through
import qchem.Symmetry.Lattice_3D.SpaceGroup;     // the library's own detection, checked against a CIF's stated ops

namespace qchem::CrystalFiles
{
namespace
{
using json = nlohmann::ordered_json;
using M3   = Matrix3D<double>;

//! Å per bohr, CODATA 2018 (the files speak Å, an entry speaks bohr).
constexpr double kAngstromPerBohr = 0.529177210903;
//! Two expanded CIF sites closer than this (Å) are ONE site.  A CIF states 4-5 decimals: a special position written 0.3333 lands
//! ~1e-3 Å from its partner after an op, while no two real atoms come within ~0.7 Å.  A tolerance between the two, not a calibration.
constexpr double kSameSiteAngstrom = 0.01;
//! SpaceGroup::Detect's position tolerance (fractional) on an imported cell: the same CIF precision, so the check is not stricter
//! than the file's own digits.
constexpr double kDetectTol = 2e-3;

struct Site { std::string el; rvec3_t f; int spin=0; };

[[noreturn]] void Refuse(const std::string& source, const std::string& what)
{
    throw std::invalid_argument(source+": "+what);
}

std::string Today()
{
    const std::time_t t=std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    std::tm tm{}; localtime_r(&t, &tm);
    std::ostringstream os; os<<std::put_time(&tm, "%Y-%m-%d");
    return os.str();
}

std::string Lower(std::string s) { for (char& c : s) c=char(std::tolower((unsigned char)c)); return s; }
std::string Trim(const std::string& s)
{
    const auto b=s.find_first_not_of(" \t\r\n"), e=s.find_last_not_of(" \t\r\n");
    return b==std::string::npos ? std::string() : s.substr(b, e-b+1);
}
std::vector<std::string> Words(const std::string& line)
{
    std::istringstream is(line); std::vector<std::string> w; std::string x;
    while (is>>x) w.push_back(x);
    return w;
}
std::vector<std::string> Lines(const std::string& text)
{
    std::vector<std::string> out; std::istringstream is(text); std::string l;
    while (std::getline(is, l)) { if (!l.empty() && l.back()=='\r') l.pop_back(); out.push_back(l); }
    return out;
}

double Number(const std::string& s, const std::string& source, const std::string& what)
{
    std::string t=s;
    if (auto p=t.find('('); p!=std::string::npos) t=t.substr(0, p);   // a CIF uncertainty: 4.1749(3)
    try { size_t used=0; const double v=std::stod(t, &used); if (used==t.size()) return v; }
    catch (const std::exception&) {}
    Refuse(source, what+" is not a number: '"+s+"'");
}

//! The element a symbol or label names: "Mg", "O1", "Na1+", "Fe_pv" -> Mg, O, Na, Fe.  Two letters first, then one.
std::string ElementOf(const std::string& token, const std::string& source)
{
    std::string letters;
    for (char c : token) { if (std::isalpha((unsigned char)c)) letters+=c; else break; }
    if (!letters.empty())
    {
        std::string two=letters.substr(0, std::min<size_t>(2, letters.size()));
        two[0]=char(std::toupper((unsigned char)two[0]));
        if (two.size()==2) two[1]=char(std::tolower((unsigned char)two[1]));
        if (two.size()==2 && thePeriodicTable().GetZ(two)>0) return two;
        const std::string one=two.substr(0, 1);
        if (thePeriodicTable().GetZ(one)>0) return one;
    }
    Refuse(source, "'"+token+"' names no element");
}

double Wrap(double x) { double w=x-std::floor(x); return w>=1.0-1e-12 ? 0.0 : w; }
rvec3_t Wrap(const rvec3_t& f) { return rvec3_t(Wrap(f.x), Wrap(f.y), Wrap(f.z)); }
//! The shortest lattice image of a fractional difference.
rvec3_t Nearest(const rvec3_t& d) { return rvec3_t(d.x-std::round(d.x), d.y-std::round(d.y), d.z-std::round(d.z)); }
double Length(const rvec3_t& v) { return std::sqrt(v.x*v.x+v.y*v.y+v.z*v.z); }

//! The entry: comment, lattice{vectors_bohr}, atoms, provenance.  A cell matrix's COLUMNS are the lattice vectors.
json MakeEntry(const M3& A, const std::vector<Site>& sites, const std::string& comment, const json& provenance)
{
    json e;
    e["comment"]=comment;
    json v=json::array();
    for (int i=1;i<=3;i++) v.push_back({A(1,i), A(2,i), A(3,i)});
    e["lattice"]={{"vectors_bohr", v}};
    json atoms=json::array();
    for (const Site& s : sites)
    {
        json a={{"el", s.el}, {"frac", {s.f.x, s.f.y, s.f.z}}};
        if (s.spin!=0) a["spin"]=s.spin;
        atoms.push_back(a);
    }
    e["atoms"]=atoms;
    e["provenance"]=provenance;
    return e;
}

std::string LegalName(const std::string& source)
{
    std::string stem=std::filesystem::path(source).stem().string();
    if (stem.empty()) stem=std::filesystem::path(source).filename().string();
    for (char& c : stem) if (!(std::isalnum((unsigned char)c) || c=='_' || c=='-')) c='_';
    if (stem.empty() || stem[0]=='_') stem="x"+stem;
    return stem;
}

std::string Formula(const std::vector<Site>& sites)
{
    std::vector<std::pair<std::string,int>> n;
    for (const Site& s : sites)
    {
        auto it=std::find_if(n.begin(), n.end(), [&](const auto& p){ return p.first==s.el; });
        if (it==n.end()) n.push_back({s.el, 1}); else it->second++;
    }
    std::string f;
    for (const auto& [el,k] : n) f+=el+(k>1 ? std::to_string(k) : std::string());
    return f;
}

//! Build the cell through THE parser (an imported entry is an entry) -- a translator that emits something it refuses is a defect.
Imported Finish(json e, const std::string& source, std::vector<std::string> warnings)
{
    if (!warnings.empty()) e["provenance"]["warnings"]=warnings;
    (void)StructureData::CellFromEntry(e.dump(), source);
    return Imported{e.dump(2), LegalName(source), std::move(warnings)};
}

json Provenance(Format f, const std::string& source, const std::string& date)
{
    return json{{"source", std::filesystem::path(source).filename().string()}, {"format", NameOf(f)},
                {"imported", date.empty() ? Today() : date}};
}

// ------------------------------------------------------------------------------------------------------------------- extended xyz

//! The comment line's key=value pairs (ASE / extxyz): values may be "double quoted"; a bare key is a flag.
std::map<std::string,std::string> KeyValues(const std::string& line)
{
    std::map<std::string,std::string> kv;
    size_t i=0;
    auto skip=[&]{ while (i<line.size() && std::isspace((unsigned char)line[i])) i++; };
    auto word=[&](bool stopAtEq)
    {
        std::string w;
        if (i<line.size() && (line[i]=='"' || line[i]=='\''))
        {
            const char q=line[i++];
            while (i<line.size() && line[i]!=q) w+=line[i++];
            if (i<line.size()) i++;
            return w;
        }
        while (i<line.size() && !std::isspace((unsigned char)line[i]) && !(stopAtEq && line[i]=='=')) w+=line[i++];
        return w;
    };
    for (skip(); i<line.size(); skip())
    {
        const std::string k=word(true);
        skip();
        if (i<line.size() && line[i]=='=') { i++; skip(); kv[Lower(k)]=word(false); }
        else if (!k.empty()) kv[Lower(k)]="T";
    }
    return kv;
}

Imported ReadExtendedXYZ(const std::string& text, const std::string& source, const std::string& date)
{
    auto L=Lines(text);
    while (!L.empty() && Trim(L.back()).empty()) L.pop_back();
    if (L.size()<2) Refuse(source, "an extended xyz needs a count line and a comment line");
    const size_t n=(size_t)Number(Trim(L[0]), source, "the atom count");
    if (n==0) Refuse(source, "the atom count is 0");
    if (L.size()<2+n) Refuse(source, "the count says "+std::to_string(n)+" atoms but the file has "+std::to_string(L.size()-2)+" atom lines");
    if (L.size()>2+n) Refuse(source, "the file holds more than one frame; an entry is one structure (split the file)");
    const auto kv=KeyValues(L[1]);
    if (!kv.count("lattice")) Refuse(source, "no Lattice=\"...\" on the comment line: a plain xyz is a MOLECULE (File > Import molecule), not a crystal");
    if (kv.count("pbc"))
    {
        const auto p=Words(kv.at("pbc"));
        const bool all=p.size()==3 && std::all_of(p.begin(), p.end(), [](const std::string& s){ const auto l=Lower(s); return l=="t" || l=="true" || l=="1"; });
        if (!all) Refuse(source, "pbc=\""+kv.at("pbc")+"\": an entry is periodic in all three directions, so a slab or a wire cannot be one");
    }
    const auto lv=Words(kv.at("lattice"));
    if (lv.size()!=9) Refuse(source, "Lattice must hold 9 numbers (a, b, c in Å), not "+std::to_string(lv.size()));
    M3 A;
    for (int i=0;i<3;i++) for (int x=0;x<3;x++) A(x+1,i+1)=Number(lv[3*i+x], source, "a Lattice component")/kAngstromPerBohr;
    if (!(std::abs(Determinant(A))>0.0)) Refuse(source, "the Lattice vectors are linearly dependent");
    // The Properties columns: name:type:count triples; species (S) and pos (R 3) are the structure, spin (I) our own export's.
    const std::string props=kv.count("properties") ? kv.at("properties") : std::string("species:S:1:pos:R:3");
    std::vector<std::string> P;
    { std::string w; std::istringstream is(props); while (std::getline(is, w, ':')) P.push_back(w); }
    if (P.size()%3!=0) Refuse(source, "Properties=\""+props+"\" is not name:type:count triples");
    int col=0, cSpecies=-1, cPos=-1, cSpin=-1;
    std::vector<std::string> ignored;
    for (size_t k=0;k<P.size();k+=3)
    {
        const std::string name=Lower(P[k]);
        const int cnt=(int)Number(P[k+2], source, "a Properties count");
        if      (name=="species") cSpecies=col;
        else if (name=="pos")     { if (cnt!=3) Refuse(source, "pos must have 3 columns"); cPos=col; }
        else if (name=="spin")    cSpin=col;
        else ignored.push_back(P[k]);
        col+=cnt;
    }
    if (cSpecies<0 || cPos<0) Refuse(source, "Properties names no species and pos columns: \""+props+"\"");
    std::vector<std::string> warnings;
    if (!ignored.empty())
    {
        std::string w="per-atom columns ignored (not structure";
        if (std::find_if(ignored.begin(), ignored.end(), [](const std::string& s){ return Lower(s).find("mag")!=std::string::npos; })!=ignored.end())
            w+="; state the site spins after import";
        w+="):";
        for (const auto& s : ignored) w+=" "+s;
        warnings.push_back(w);
    }
    const M3 Ainv=Invert(A);
    std::vector<Site> sites;
    for (size_t i=0;i<n;i++)
    {
        const auto w=Words(L[2+i]);
        if ((int)w.size()<col) Refuse(source, "atom line "+std::to_string(i+1)+" has "+std::to_string(w.size())+" columns, Properties says "+std::to_string(col));
        const rvec3_t r(Number(w[cPos], source, "x")/kAngstromPerBohr, Number(w[cPos+1], source, "y")/kAngstromPerBohr,
                        Number(w[cPos+2], source, "z")/kAngstromPerBohr);
        Site s{ElementOf(w[cSpecies], source), Ainv*r};
        if (cSpin>=0)
        {
            s.spin=(int)Number(w[cSpin], source, "a spin");
            if (s.spin<-1 || s.spin>1) Refuse(source, "a site spin is +1, -1 or 0, not "+w[cSpin]);
        }
        sites.push_back(s);
    }
    json prov=Provenance(Format::ExtendedXYZ, source, date);
    const std::string comment="imported from "+prov["source"].get<std::string>()+" (extended xyz): "+Formula(sites)
                             +(kv.count("comment") ? ", "+kv.at("comment") : std::string());
    return Finish(MakeEntry(A, sites, comment, prov), source, warnings);
}

// ------------------------------------------------------------------------------------------------------------------------- POSCAR

Imported ReadPOSCAR(const std::string& text, const std::string& source, const std::string& date)
{
    const auto L=Lines(text);
    if (L.size()<8) Refuse(source, "a POSCAR needs at least 8 lines (comment, scale, 3 vectors, species, counts, mode)");
    const auto sc=Words(L[1]);
    if (sc.size()!=1 && sc.size()!=3) Refuse(source, "line 2 (the scale) must hold 1 or 3 numbers");
    M3 A;   // Å until scaled
    for (int i=0;i<3;i++)
    {
        const auto w=Words(L[2+i]);
        if (w.size()<3) Refuse(source, "lattice line "+std::to_string(i+1)+" needs 3 numbers");
        for (int x=0;x<3;x++) A(x+1,i+1)=Number(w[x], source, "a lattice component");
    }
    rvec3_t s3(1,1,1);
    if (sc.size()==3) s3=rvec3_t(Number(sc[0], source, "a scale"), Number(sc[1], source, "a scale"), Number(sc[2], source, "a scale"));
    else
    {
        const double s=Number(sc[0], source, "the scale");
        if (s<0.0) { const double v=std::abs(Determinant(A)); const double k=std::cbrt(-s/v); s3=rvec3_t(k,k,k); }   // negative = the cell volume, Å^3
        else s3=rvec3_t(s,s,s);
    }
    for (int i=1;i<=3;i++) for (int x=1;x<=3;x++) A(x,i)*= (x==1 ? s3.x : x==2 ? s3.y : s3.z) / kAngstromPerBohr;
    if (!(std::abs(Determinant(A))>0.0)) Refuse(source, "the lattice vectors are linearly dependent");
    size_t line=5;
    auto names=Words(L[line]);
    if (names.empty() || !std::isalpha((unsigned char)names[0][0]))
        Refuse(source, "line 6 holds counts, not element names: a VASP 4 POSCAR names no species -- add the element names above the counts");
    for (auto& nm : names) nm=ElementOf(nm.substr(0, nm.find_first_of("_/")), source);
    const auto counts=Words(L[++line]);
    if (counts.size()!=names.size()) Refuse(source, std::to_string(names.size())+" species names but "+std::to_string(counts.size())+" counts");
    std::vector<std::string> el;
    for (size_t k=0;k<counts.size();k++) for (int j=0, c=(int)Number(counts[k], source, "a count"); j<c; j++) el.push_back(names[k]);
    std::vector<std::string> warnings;
    std::string mode=Trim(L[++line]);
    if (!mode.empty() && (mode[0]=='s' || mode[0]=='S'))
    {
        warnings.push_back("selective-dynamics flags ignored (a relaxation input, not structure)");
        if (++line>=L.size()) Refuse(source, "the file ends after 'Selective dynamics'");
        mode=Trim(L[line]);
    }
    const char m=mode.empty() ? ' ' : char(std::tolower((unsigned char)mode[0]));
    if (m!='d' && m!='c' && m!='k') Refuse(source, "line "+std::to_string(line+1)+" must say Direct or Cartesian, not '"+mode+"'");
    if (L.size()<line+1+el.size()) Refuse(source, "the counts say "+std::to_string(el.size())+" atoms but fewer position lines follow");
    const M3 Ainv=Invert(A);
    std::vector<Site> sites;
    for (size_t i=0;i<el.size();i++)
    {
        const auto w=Words(L[line+1+i]);
        if (w.size()<3) Refuse(source, "position line "+std::to_string(i+1)+" needs 3 numbers");
        const rvec3_t p(Number(w[0], source, "a coordinate"), Number(w[1], source, "a coordinate"), Number(w[2], source, "a coordinate"));
        // Cartesian positions are in Å, scaled by the same factor as the lattice
        const rvec3_t f= m=='d' ? p : Ainv*rvec3_t(p.x*s3.x/kAngstromPerBohr, p.y*s3.y/kAngstromPerBohr, p.z*s3.z/kAngstromPerBohr);
        sites.push_back(Site{el[i], f});
    }
    json prov=Provenance(Format::POSCAR, source, date);
    const std::string title=Trim(L[0]);
    const std::string comment="imported from "+prov["source"].get<std::string>()+" (POSCAR): "+Formula(sites)+(title.empty() ? "" : ", "+title);
    return Finish(MakeEntry(A, sites, comment, prov), source, warnings);
}

// ---------------------------------------------------------------------------------------------------------------------------- CIF

//! One CIF data block: its single-valued tags and its loops (tags lower-cased: CIF tags are case-insensitive).
struct CifBlock
{
    std::string name;
    std::map<std::string,std::string> tags;
    std::vector<std::pair<std::vector<std::string>, std::vector<std::vector<std::string>>>> loops;   // (tags, rows)

    std::optional<std::string> Tag(std::initializer_list<const char*> keys) const
    {
        for (const char* k : keys) if (auto it=tags.find(k); it!=tags.end() && it->second!="?" && it->second!=".") return it->second;
        return std::nullopt;
    }
    //! The loop holding \a tag, and that tag's column.
    std::optional<std::pair<size_t,size_t>> Column(const std::string& tag) const
    {
        for (size_t l=0;l<loops.size();l++)
            for (size_t c=0;c<loops[l].first.size();c++) if (loops[l].first[c]==tag) return std::make_pair(l, c);
        return std::nullopt;
    }
};

//! The CIF 1.1 token stream: whitespace-separated words, 'single' / "double" quoted strings (a quote closes only before
//! whitespace), ;-delimited text fields, # comments.
std::vector<std::string> CifTokens(const std::string& text)
{
    std::vector<std::string> tok;
    const auto L=Lines(text);
    for (size_t li=0; li<L.size(); li++)
    {
        const std::string& l=L[li];
        if (!l.empty() && l[0]==';')
        {   // a text field: up to the next line that starts with ';'
            std::string field=l.substr(1);
            while (++li<L.size() && !(!L[li].empty() && L[li][0]==';')) field+="\n"+L[li];
            tok.push_back(Trim(field));
            continue;
        }
        size_t i=0;
        while (i<l.size())
        {
            while (i<l.size() && std::isspace((unsigned char)l[i])) i++;
            if (i>=l.size() || l[i]=='#') break;
            if (l[i]=='\'' || l[i]=='"')
            {
                const char q=l[i++];
                std::string w;
                while (i<l.size() && !(l[i]==q && (i+1==l.size() || std::isspace((unsigned char)l[i+1])))) w+=l[i++];
                i++;
                tok.push_back(w);
                continue;
            }
            std::string w;
            while (i<l.size() && !std::isspace((unsigned char)l[i])) w+=l[i++];
            tok.push_back(w);
        }
    }
    return tok;
}

std::vector<CifBlock> ParseCif(const std::string& text, const std::string& source)
{
    const auto T=CifTokens(text);
    std::vector<CifBlock> blocks;
    auto isTag=[](const std::string& t){ return !t.empty() && t[0]=='_'; };
    auto isKeyword=[](const std::string& t){ const auto l=Lower(t); return l=="loop_" || l.rfind("data_",0)==0 || l.rfind("save_",0)==0 || l=="global_"; };
    for (size_t i=0; i<T.size();)
    {
        const std::string low=Lower(T[i]);
        if (low.rfind("data_",0)==0) { blocks.push_back(CifBlock{T[i].substr(5)}); i++; continue; }
        if (blocks.empty()) Refuse(source, "CIF content before any data_ block");
        CifBlock& b=blocks.back();
        if (low=="loop_")
        {
            std::vector<std::string> tags;
            for (i++; i<T.size() && isTag(T[i]); i++) tags.push_back(Lower(T[i]));
            std::vector<std::string> vals;
            for (; i<T.size() && !isTag(T[i]) && !isKeyword(T[i]); i++) vals.push_back(T[i]);
            if (tags.empty()) Refuse(source, "a loop_ with no tags");
            if (vals.size()%tags.size()!=0)
                Refuse(source, "a loop over "+tags[0]+"... holds "+std::to_string(vals.size())+" values, not a multiple of its "+std::to_string(tags.size())+" tags");
            std::vector<std::vector<std::string>> rows;
            for (size_t k=0;k<vals.size();k+=tags.size()) rows.emplace_back(vals.begin()+k, vals.begin()+k+tags.size());
            b.loops.push_back({tags, rows});
            continue;
        }
        if (isTag(T[i]))
        {
            if (i+1>=T.size() || isTag(T[i+1]) || isKeyword(T[i+1])) Refuse(source, "tag "+T[i]+" has no value");
            b.tags[Lower(T[i])]=T[i+1];
            i+=2;
            continue;
        }
        i++;   // save_ frames, stray words: not structure
    }
    return blocks;
}

struct SymOp { M3 W; rvec3_t tau; };

//! "x,1/2+y,-z+1/2" / "-x+y, -x, z+1/3" / "+x,+y,+z" -> {W|tau} in the CIF cell's fractional basis.
SymOp ParseSymOp(const std::string& s, const std::string& source)
{
    std::vector<std::string> parts;
    { std::string p; std::istringstream is(s); while (std::getline(is, p, ',')) parts.push_back(p); }
    if (parts.size()!=3) Refuse(source, "symmetry operation '"+s+"' does not have three components");
    SymOp op; op.W=M3(0,0,0, 0,0,0, 0,0,0);
    double tau[3]={0,0,0};
    for (int r=0;r<3;r++)
    {
        std::string e; for (char c : parts[r]) if (!std::isspace((unsigned char)c)) e+=char(std::tolower((unsigned char)c));
        size_t i=0;
        while (i<e.size())
        {
            double sign=1.0;
            while (i<e.size() && (e[i]=='+' || e[i]=='-')) { if (e[i]=='-') sign=-sign; i++; }
            std::string num;
            while (i<e.size() && (std::isdigit((unsigned char)e[i]) || e[i]=='.' || e[i]=='/')) num+=e[i++];
            if (i<e.size() && e[i]=='*') i++;
            double v=1.0;
            if (!num.empty())
            {
                const auto sl=num.find('/');
                v= sl==std::string::npos ? Number(num, source, "a symop number")
                                         : Number(num.substr(0, sl), source, "a symop fraction")/Number(num.substr(sl+1), source, "a symop fraction");
            }
            if (i<e.size() && (e[i]=='x' || e[i]=='y' || e[i]=='z')) { op.W(r+1, 1+(e[i]-'x'))+=sign*v; i++; }
            else if (!num.empty()) tau[r]+=sign*v;
            else Refuse(source, "symmetry operation '"+s+"': cannot read '"+e+"'");
        }
    }
    op.tau=rvec3_t(tau[0], tau[1], tau[2]);
    return op;
}

bool IsIdentityW(const M3& W)
{
    for (int i=1;i<=3;i++) for (int j=1;j<=3;j++) if (std::abs(W(i,j)-(i==j ? 1.0 : 0.0))>1e-9) return false;
    return true;
}

using IntW = std::array<int,9>;
std::optional<IntW> Integer(const M3& W)
{
    IntW k;
    for (int i=1;i<=3;i++) for (int j=1;j<=3;j++)
    {
        const double r=std::round(W(i,j));
        if (std::abs(W(i,j)-r)>1e-6) return std::nullopt;
        k[3*(i-1)+j-1]=(int)r;
    }
    return k;
}

//! THE PRIMITIVE CELL of a centred CIF lattice: the lattice generated by the unit vectors and the centring translations
//! \a centring (fractional, of the conventional cell) has a basis of three of its shortest vectors with |det| = 1/n
//! (n = 1 + the number of centrings).  Returns P, its columns those vectors (f_conv = P f_prim).
M3 PrimitiveBasis(const M3& A, const std::vector<rvec3_t>& centring, const std::string& source)
{
    const double n=1.0+double(centring.size());
    std::vector<rvec3_t> cand;
    std::vector<rvec3_t> shifts={rvec3_t(0,0,0)}; for (const auto& t : centring) shifts.push_back(t);
    for (const auto& t : shifts)
        for (int a=-1;a<=1;a++) for (int b=-1;b<=1;b++) for (int c=-1;c<=1;c++)
        {
            const rvec3_t v(t.x+a, t.y+b, t.z+c);
            if (Length(v)>1e-9) cand.push_back(v);
        }
    std::sort(cand.begin(), cand.end(), [&](const rvec3_t& u, const rvec3_t& v){ return Length(A*u)<Length(A*v); });
    // Among the triples that span the lattice, the SHORTEST, and of equal length the one with the fewest negative components:
    // the textbook choice ((0,1/2,1/2), (1/2,0,1/2), (1/2,1/2,0) for F; 60 degrees, not 120), so an import is reproducible.
    const size_t N=std::min<size_t>(cand.size(), 60);
    auto negatives=[](const rvec3_t& v){ return int(v.x<-1e-9)+int(v.y<-1e-9)+int(v.z<-1e-9); };
    std::optional<M3> best; double bestLen=0.0; int bestNeg=0;
    for (size_t i=0;i<N;i++) for (size_t j=i+1;j<N;j++) for (size_t k=j+1;k<N;k++)
    {
        M3 P(cand[i].x, cand[j].x, cand[k].x, cand[i].y, cand[j].y, cand[k].y, cand[i].z, cand[j].z, cand[k].z);
        const double d=Determinant(P);
        if (std::abs(std::abs(d)-1.0/n)>1e-6) continue;
        const double len=Length(A*cand[i])+Length(A*cand[j])+Length(A*cand[k]);
        const int neg=negatives(cand[i])+negatives(cand[j])+negatives(cand[k]);
        if (best && (len>bestLen*(1.0+1e-9) || (len>bestLen*(1.0-1e-9) && neg>=bestNeg))) continue;
        if (d<0.0) { P=M3(cand[i].x, cand[k].x, cand[j].x, cand[i].y, cand[k].y, cand[j].y, cand[i].z, cand[k].z, cand[j].z); }   // right-handed
        best=P; bestLen=len; bestNeg=neg;
    }
    if (best) return *best;
    Refuse(source, "its "+std::to_string(centring.size())+" centring translation(s) generate no primitive cell (a CIF symop list that is not a group?)");
}

std::string Stated(const CifBlock& b)
{
    return Trim(b.Tag({"_space_group_name_h-m_alt", "_symmetry_space_group_name_h-m", "_space_group_name_hall", "_symmetry_space_group_name_hall"}).value_or(""));
}

Imported ReadCIF(const std::string& text, const std::string& source, const std::string& date)
{
    auto blocks=ParseCif(text, source);
    std::vector<CifBlock> withAtoms;
    for (auto& b : blocks) if (b.Column("_atom_site_fract_x")) withAtoms.push_back(b);
    if (withAtoms.empty())
    {
        for (auto& b : blocks) if (b.Column("_atom_site_cartn_x")) Refuse(source, "its sites are Cartesian (_atom_site_Cartn_x) only; an import reads fractional sites");
        Refuse(source, "no data block lists fractional sites (_atom_site_fract_x)");
    }
    if (withAtoms.size()>1)
    {
        std::string names; for (const auto& b : withAtoms) names+=(names.empty() ? "" : ", ")+b.name;
        Refuse(source, "it holds "+std::to_string(withAtoms.size())+" structures (data blocks "+names+"); an entry is one -- split the file");
    }
    const CifBlock& b=withAtoms[0];
    auto cell=[&](const char* tag)
    {
        const auto v=b.Tag({tag});
        if (!v) Refuse(source, std::string("no ")+tag);
        return Number(*v, source, tag);
    };
    const double a=cell("_cell_length_a"), bb=cell("_cell_length_b"), c=cell("_cell_length_c");
    const double al=cell("_cell_angle_alpha"), be=cell("_cell_angle_beta"), ga=cell("_cell_angle_gamma");
    const M3 A=CellGeometry(a/kAngstromPerBohr, bb/kAngstromPerBohr, c/kAngstromPerBohr, al, be, ga).GetCellMatrix();

    // THE SYMMETRY OPERATIONS the file states.  A file naming a group without listing its ops cannot be expanded: the library
    // carries no space-group tables, and guessing a setting from a symbol is exactly the silent substitution to avoid.
    std::vector<std::string> warnings;
    std::vector<SymOp> ops;
    const std::string statedName=Stated(b);
    std::optional<std::pair<size_t,size_t>> sc=b.Column("_space_group_symop_operation_xyz");
    if (!sc) sc=b.Column("_symmetry_equiv_pos_as_xyz");
    if (sc) for (const auto& row : b.loops[sc->first].second) ops.push_back(ParseSymOp(row[sc->second], source));
    else if (auto one=b.Tag({"_space_group_symop_operation_xyz", "_symmetry_equiv_pos_as_xyz"})) ops.push_back(ParseSymOp(*one, source));
    if (ops.empty())
    {
        if (!statedName.empty() && Lower(statedName)!="p 1" && Lower(statedName)!="p1")
            Refuse(source, "it names space group '"+statedName+"' but lists no symmetry operations; the library carries no space-group tables to "
                           "expand it -- take the CIF with its symop loop (the database's own download has it)");
        warnings.push_back("no symmetry operations listed: read as P 1 (the sites as given)");
        ops.push_back(SymOp{M3(), rvec3_t(0,0,0)});
    }
    if (std::none_of(ops.begin(), ops.end(), [](const SymOp& o){ return IsIdentityW(o.W) && Length(Nearest(o.tau))<1e-9; }))
        Refuse(source, "its symmetry operations do not include x,y,z");

    // THE ASYMMETRIC UNIT.  What an ordered entry cannot represent is refused by name: a partial occupancy, a disorder group.
    const auto [lAtoms, cx]=*b.Column("_atom_site_fract_x");
    const auto& loop=b.loops[lAtoms];
    auto col=[&](const char* tag)->std::optional<size_t>
    { for (size_t k=0;k<loop.first.size();k++) if (loop.first[k]==tag) return k; return std::nullopt; };
    const auto cy=col("_atom_site_fract_y"), cz=col("_atom_site_fract_z"), cLabel=col("_atom_site_label"), cType=col("_atom_site_type_symbol");
    const auto cOcc=col("_atom_site_occupancy"), cDis=col("_atom_site_disorder_group"), cAss=col("_atom_site_disorder_assembly");
    const auto cCalc=col("_atom_site_calc_flag");
    if (!cy || !cz) Refuse(source, "the atom loop lacks _atom_site_fract_y / _atom_site_fract_z");
    if (!cLabel && !cType) Refuse(source, "the atom loop names no element (no _atom_site_type_symbol or _atom_site_label)");
    struct Asym { std::string label, el; rvec3_t f; };
    std::vector<Asym> asym;
    for (const auto& row : loop.second)
    {
        const std::string label=cLabel ? row[*cLabel] : row[*cType];
        if (cCalc && Lower(row[*cCalc])=="dum") continue;   // a dummy site (a centroid), not an atom
        if (cOcc && row[*cOcc]!="?" && row[*cOcc]!=".")
        {
            const double occ=Number(row[*cOcc], source, "site "+label+"'s occupancy");
            if (std::abs(occ-1.0)>1e-3)
                Refuse(source, "site "+label+" has occupancy "+row[*cOcc]+": a materials entry is an ORDERED crystal and cannot hold a partially "
                               "occupied site (build an ordered supercell for the composition and import that)");
        }
        for (auto cd : {cDis, cAss})
            if (cd && row[*cd]!="." && row[*cd]!="?")
                Refuse(source, "site "+label+" belongs to disorder group / assembly '"+row[*cd]+"': a materials entry is an ORDERED crystal and cannot hold "
                               "alternative positions (choose one ordering and import that)");
        const std::string el=ElementOf(cType && row[*cType]!="?" && row[*cType]!="." ? row[*cType] : label, source);
        asym.push_back({label, el, Wrap(rvec3_t(Number(row[cx], source, label+" x"), Number(row[*cy], source, label+" y"), Number(row[*cz], source, label+" z")))});
    }
    if (asym.empty()) Refuse(source, "the atom loop lists no atoms");

    // SYMMETRY-EXPAND the asymmetric unit onto the conventional cell.  Two images within kSameSiteAngstrom are one site; two
    // DIFFERENT sites landing there is refused (it is disorder by another name, or a wrong setting).
    const double same=kSameSiteAngstrom/kAngstromPerBohr;
    std::vector<Site> conv;
    std::vector<std::string> convLabel;
    for (const Asym& s : asym)
        for (const SymOp& o : ops)
        {
            const rvec3_t f=Wrap(o.W*s.f+o.tau);
            bool dup=false;
            for (size_t k=0;k<conv.size() && !dup;k++)
                if (Length(A*Nearest(f-conv[k].f))<same)
                {
                    if (convLabel[k]!=s.label) Refuse(source, "sites "+convLabel[k]+" and "+s.label+" fall on one point after the symmetry expansion");
                    dup=true;
                }
            if (!dup) { conv.push_back(Site{s.el, f}); convLabel.push_back(s.label); }
        }

    // THE PRIMITIVE CELL of a centred lattice, from the centring translations the file itself lists.
    std::vector<rvec3_t> centring;
    for (const SymOp& o : ops)
        if (IsIdentityW(o.W) && Length(Nearest(o.tau))>1e-6)
        {
            const rvec3_t t=Wrap(o.tau);
            if (std::none_of(centring.begin(), centring.end(), [&](const rvec3_t& u){ return Length(Nearest(u-t))<1e-6; })) centring.push_back(t);
        }
    M3 P;   // identity
    M3 Ap=A;
    std::vector<Site> sites=conv;
    std::string cellNote="the CIF's own cell";
    if (!centring.empty())
    {
        P=PrimitiveBasis(A, centring, source);
        Ap=A*P;
        const M3 Pinv=Invert(P);
        sites.clear();
        for (const Site& s : conv)
        {
            const rvec3_t f=Wrap(Pinv*s.f);
            if (std::none_of(sites.begin(), sites.end(), [&](const Site& u){ return Length(Ap*Nearest(u.f-f))<same; })) sites.push_back(Site{s.el, f});
        }
        if (sites.size()*(centring.size()+1)!=conv.size())
            Refuse(source, "the "+std::to_string(conv.size())+"-site conventional cell does not reduce to a primitive one under its "
                           +std::to_string(centring.size())+" centring translation(s) (the expanded sites are not centring-invariant)");
        cellNote="the primitive cell of the CIF's "+std::to_string(conv.size())+"-site conventional cell ("+std::to_string(centring.size()+1)+" lattice points in it)";
    }

    // THE CHECK: the stated ops' linear parts (in the stored cell's basis) against the library's own detection on the stored cell.
    std::set<IntW> statedW;
    bool latticeCompatible=true;
    const M3 Pinv=Invert(P);
    for (const SymOp& o : ops)
    {
        if (auto k=Integer(Pinv*o.W*P)) statedW.insert(*k); else latticeCompatible=false;
    }
    std::vector<Symmetry::Lattice_3D::AtomSite> basis;
    for (const Site& s : sites) basis.push_back({(int)thePeriodicTable().GetZ(s.el), s.f});
    const auto detected=Symmetry::Lattice_3D::SpaceGroup::Detect(Ap, basis, kDetectTol);
    std::set<IntW> detectedW;
    for (const auto& op : detected.Ops()) if (auto k=Integer(op.W)) detectedW.insert(*k);
    const bool statedInDetected=std::includes(detectedW.begin(), detectedW.end(), statedW.begin(), statedW.end());
    const bool agrees=latticeCompatible && statedInDetected && statedW.size()==detectedW.size();
    json sg={{"stated", statedName}, {"statedPointOps", (long)statedW.size()}, {"detectedPointOps", (long)detectedW.size()}, {"agrees", agrees}};
    if (auto num=b.Tag({"_space_group_it_number", "_symmetry_int_tables_number"})) sg["number"]=(long)Number(*num, source, "the space-group number");
    const std::string sgName=statedName.empty() ? std::string("the stated operations") : "'"+statedName+"'";
    if (!latticeCompatible || !statedInDetected)
        warnings.push_back("space group: "+sgName+" lists operations the library does NOT find in the expanded structure ("+std::to_string(statedW.size())
                           +" stated, "+std::to_string(detectedW.size())+" detected point operations): check the coordinates' precision or the setting");
    else if (detectedW.size()>statedW.size() && statedW.size()>1)   // a P 1 file states no group: nothing to under-state
        warnings.push_back("space group: the library detects MORE symmetry than "+sgName+" states ("+std::to_string(detectedW.size())+" point operations against "
                           +std::to_string(statedW.size())+"): the file under-states its group, or a distortion is below the coordinates' precision");

    json prov=Provenance(Format::CIF, source, date);
    if (auto cod=b.Tag({"_cod_database_code"}))      prov["database"]="COD "+*cod;
    else if (auto ic=b.Tag({"_database_code_icsd"})) prov["database"]="ICSD "+*ic;
    else if (auto am=b.Tag({"_database_code_amcsd"})) prov["database"]="AMCSD "+*am;
    if (auto doi=b.Tag({"_journal_paper_doi", "_citation_doi"})) prov["doi"]=*doi;
    else if (auto dc=b.Column("_citation_doi"); dc && !b.loops[dc->first].second.empty()) prov["doi"]=b.loops[dc->first].second[0][dc->second];
    prov["dataBlock"]=b.name;
    prov["cell"]=cellNote;
    prov["spaceGroup"]=sg;
    const std::string what=b.Tag({"_chemical_name_mineral", "_chemical_name_common", "_chemical_name_systematic", "_chemical_name"}).value_or("");
    std::ostringstream cm;
    cm<<"imported from "<<prov["source"].get<std::string>()<<" (CIF"<<(prov.contains("database") ? ", "+prov["database"].get<std::string>() : std::string())<<"): "
      <<(what.empty() ? std::string() : what+", ")<<Formula(sites)<<(statedName.empty() ? std::string() : ", "+statedName)
      <<", a = "<<a<<" Å";
    return Finish(MakeEntry(Ap, sites, cm.str(), prov), source, warnings);
}

// ------------------------------------------------------------------------------------------------------------------------ writers

std::string Fixed(double x) { std::ostringstream os; os<<std::setprecision(17)<<x; return os.str(); }

} // anonymous

std::string NameOf(Format f)
{
    switch (f)
    {
        case Format::CIF:         return "CIF";
        case Format::POSCAR:      return "POSCAR";
        case Format::ExtendedXYZ: return "extended xyz";
    }
    return "?";
}

Format FormatOf(const std::string& fileName)
{
    const std::filesystem::path p(fileName);
    const std::string ext=Lower(p.extension().string()), base=Lower(p.filename().string());
    if (ext==".cif") return Format::CIF;
    if (ext==".xyz" || ext==".extxyz") return Format::ExtendedXYZ;
    if (ext==".vasp" || ext==".poscar" || base.rfind("poscar",0)==0 || base.rfind("contcar",0)==0) return Format::POSCAR;
    throw std::invalid_argument(fileName+": not a crystal file the import reads (.cif; .xyz / .extxyz; POSCAR / CONTCAR / .vasp / .poscar)");
}

Imported Read(const std::string& text, Format f, const std::string& source, const std::string& date)
{
    switch (f)
    {
        case Format::CIF:         return ReadCIF(text, source, date);
        case Format::POSCAR:      return ReadPOSCAR(text, source, date);
        case Format::ExtendedXYZ: return ReadExtendedXYZ(text, source, date);
    }
    throw std::logic_error("CrystalFiles::Read: unknown format");
}

Imported ReadFile(const std::string& path, const std::string& date)
{
    const Format f=FormatOf(path);
    std::ifstream in(path, std::ios::binary);
    if (!in) throw std::invalid_argument(path+": cannot open");
    std::ostringstream os; os<<in.rdbuf();
    return Read(os.str(), f, path, date);
}

std::string Write(const std::string& entryJson, Format f, const std::string& name)
{
    json e;
    try { e=json::parse(entryJson); }
    catch (const json::exception& x) { throw std::invalid_argument("CrystalFiles::Write: '"+name+"' is not JSON: "+x.what()); }
    const M3 A=StructureData::CellFromEntry(entryJson, name).GetCellMatrix();   // a Bravais entry's vectors, from the one parser
    std::vector<Site> sites;
    for (const auto& a : e.at("atoms"))
    {
        const auto& fr=a.at("frac");
        sites.push_back(Site{a.at("el").get<std::string>(), rvec3_t(fr[0].get<double>(), fr[1].get<double>(), fr[2].get<double>()),
                             a.contains("spin") ? a["spin"].get<int>() : 0});
    }
    const bool spins=std::any_of(sites.begin(), sites.end(), [](const Site& s){ return s.spin!=0; });
    std::string spinNote;
    if (spins) { spinNote="site spins (site order):"; for (const Site& s : sites) spinNote+=" "+std::to_string(s.spin); }
    std::ostringstream os;
    auto ang=[&](double bohr){ return Fixed(bohr*kAngstromPerBohr); };
    switch (f)
    {
    case Format::ExtendedXYZ:
    {
        os<<sites.size()<<"\n"<<"Lattice=\"";
        for (int i=1;i<=3;i++) for (int x=1;x<=3;x++) os<<(i+x>2 ? " " : "")<<ang(A(x,i));
        os<<"\" Properties=species:S:1:pos:R:3"<<(spins ? ":spin:I:1" : "")<<" pbc=\"T T T\" comment=\""<<name<<"\"\n";
        for (const Site& s : sites)
        {
            const rvec3_t r=A*s.f;
            os<<s.el<<" "<<ang(r.x)<<" "<<ang(r.y)<<" "<<ang(r.z);
            if (spins) os<<" "<<s.spin;
            os<<"\n";
        }
        break;
    }
    case Format::POSCAR:
    {
        os<<name<<(spins ? "  # "+spinNote : std::string())<<"\n1.0\n";
        for (int i=1;i<=3;i++) os<<"  "<<ang(A(1,i))<<" "<<ang(A(2,i))<<" "<<ang(A(3,i))<<"\n";
        // species in RUNS, in site order (a repeated name is a new run): the file's order is the entry's, so the round trip is exact
        std::vector<std::pair<std::string,int>> runs;
        for (const Site& s : sites) { if (runs.empty() || runs.back().first!=s.el) runs.push_back({s.el, 0}); runs.back().second++; }
        for (const auto& r : runs) os<<" "<<r.first;
        os<<"\n";
        for (const auto& r : runs) os<<" "<<r.second;
        os<<"\nDirect\n";
        for (const Site& s : sites) os<<"  "<<Fixed(s.f.x)<<" "<<Fixed(s.f.y)<<" "<<Fixed(s.f.z)<<"\n";
        break;
    }
    case Format::CIF:
    {
        // The cell parameters and EVERY site, in P 1: an entry states no group, and listing the detected operations would make the
        // import re-expand each site onto the others (a new site order).  A P 1 file asserts nothing, so its import checks nothing.
        auto len=[&](int i){ return std::sqrt(A(1,i)*A(1,i)+A(2,i)*A(2,i)+A(3,i)*A(3,i)); };
        auto angle=[&](int i, int j){ double d=0; for (int x=1;x<=3;x++) d+=A(x,i)*A(x,j); return std::acos(d/(len(i)*len(j)))*180.0/3.14159265358979323846; };
        std::string block=name; for (char& c : block) if (std::isspace((unsigned char)c)) c='_';
        os<<"data_"<<block<<"\n";
        if (e.contains("comment")) os<<"# "<<e["comment"].get<std::string>()<<"\n";
        if (spins) os<<"# "<<spinNote<<"\n";
        os<<"_cell_length_a "<<ang(len(1))<<"\n_cell_length_b "<<ang(len(2))<<"\n_cell_length_c "<<ang(len(3))<<"\n"
          <<"_cell_angle_alpha "<<Fixed(angle(2,3))<<"\n_cell_angle_beta "<<Fixed(angle(1,3))<<"\n_cell_angle_gamma "<<Fixed(angle(1,2))<<"\n"
          <<"_space_group_name_H-M_alt 'P 1'\nloop_\n_space_group_symop_operation_xyz\n'x,y,z'\n";
        os<<"loop_\n_atom_site_label\n_atom_site_type_symbol\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n_atom_site_occupancy\n";
        for (size_t k=0;k<sites.size();k++)
            os<<sites[k].el<<k+1<<" "<<sites[k].el<<" "<<Fixed(sites[k].f.x)<<" "<<Fixed(sites[k].f.y)<<" "<<Fixed(sites[k].f.z)<<" 1\n";
        break;
    }
    }
    return os.str();
}

} // namespace qchem::CrystalFiles
