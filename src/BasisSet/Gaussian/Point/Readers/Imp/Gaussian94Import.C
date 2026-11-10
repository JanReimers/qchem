// File: BasisSet/Gaussian/Point/Readers/Imp/Gaussian94Import.C  The translator itself (see the interface for what it rewrites).
module;
#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <set>
#include <sstream>
#include <string>
#include <vector>
module qchem.BasisSet.Gaussian.Point.Readers.Gaussian94Import;
import qchem.PeriodicTable;

namespace qchem::BasisSet::Gaussian
{
namespace
{
constexpr const char* kLabels = "SPDFGHI";   //!< the angular labels the Gaussian94Reader builds (l = 0..6)

std::vector<std::string> Tokens(const std::string& line)
{
    std::istringstream is(line);
    std::vector<std::string> t;
    for (std::string w; is >> w;) t.push_back(w);
    return t;
}

//! A Fortran real ("1.30D+01", "0.5D-3") as C text ("1.30E+01"); false when it is not a number.
bool CNumber(std::string s, std::string& out)
{
    for (char& c : s) if (c=='D' || c=='d') c='E';
    char* end=nullptr;
    (void)std::strtod(s.c_str(), &end);
    if (s.empty() || end!=s.c_str()+s.size()) return false;
    out=s;
    return true;
}

//! The provenance preamble may hold anything except the reader's "!" marker token: every '!' goes.
std::string Preamble(std::string s)
{
    s.erase(std::remove(s.begin(), s.end(), '!'), s.end());
    while (!s.empty() && std::isspace((unsigned char)s.front())) s.erase(s.begin());
    while (!s.empty() && std::isspace((unsigned char)s.back())) s.pop_back();
    return s;
}

struct Shell
{
    std::string                           labels;    //!< "S", "SP", "D", ...
    std::vector<std::string>              exponents;
    std::vector<std::vector<std::string>> coeffs;    //!< [primitive][label]
};
}

Outcome<ImportedBasis, std::string> ImportGaussian94(const std::string& text, const std::string& title,
                                                     const std::vector<std::string>& provenance)
{
    using O=Outcome<ImportedBasis, std::string>;
    std::vector<std::string> lines;
    {
        std::istringstream is(text);
        for (std::string l; std::getline(is, l);)
        {
            if (!l.empty() && l.back()=='\r') l.pop_back();
            lines.push_back(l);
        }
    }

    ImportedBasis out;
    std::vector<std::string> comments;
    std::ostringstream body;
    std::set<std::string> seen;
    bool hadD=false;
    int splitSP=0;
    size_t i=0;
    auto blank=[](const std::string& l){ return std::all_of(l.begin(), l.end(), [](unsigned char c){ return std::isspace(c); }); };
    auto where=[&](size_t k){ return " (line "+std::to_string(k+1)+")"; };

    while (i<lines.size())
    {
        const std::string& l=lines[i];
        if (blank(l)) { ++i; continue; }
        const auto t=Tokens(l);
        if (t[0][0]=='!') { comments.push_back(Preamble(l)); ++i; continue; }
        if (l.find("-ECP")!=std::string::npos || t[0]=="ECP")
            return O::Fail("an ECP block"+where(i)+": the library has no route for effective core potentials yet -- choose an "
                           "all-electron set");
        // An element block: "<symbol> 0".
        if (t.size()!=2 || t[1]!="0")
            return O::Fail("expected an element header '<symbol> 0'"+where(i)+", found '"+l+"'");
        std::string sym=t[0];
        for (char& c : sym) c=(char)std::tolower((unsigned char)c);
        sym[0]=(char)std::toupper((unsigned char)sym[0]);
        if (thePeriodicTable().GetZ(sym)==0) return O::Fail("'"+t[0]+"' is not an element"+where(i));
        std::string SYM=sym;
        for (char& c : SYM) c=(char)std::toupper((unsigned char)c);
        if (!seen.insert(SYM).second) return O::Fail("element "+sym+" appears twice"+where(i));
        out.elements.push_back(SYM);
        ++i;

        std::vector<Shell> shells;
        for (;;)
        {
            if (i>=lines.size()) return O::Fail("element "+sym+" has no '****' terminator");
            if (blank(lines[i])) { ++i; continue; }
            const auto h=Tokens(lines[i]);
            if (h[0]=="****") { ++i; break; }
            if (lines[i].find("-ECP")!=std::string::npos)
                return O::Fail("an ECP block in element "+sym+where(i)+": the library has no route for effective core potentials yet");
            if (h.size()!=3) return O::Fail("element "+sym+": expected a shell header '<L> <n> <scale>'"+where(i)+", found '"+lines[i]+"'");
            Shell s;
            for (char c : h[0])
            {
                const char* p=std::strchr(kLabels, std::toupper((unsigned char)c));
                if (!p || !*p) return O::Fail("element "+sym+": angular label '"+h[0]+"'"+where(i)+" -- the reader builds S through I");
                s.labels+=(char)std::toupper((unsigned char)c);
                out.maxL=std::max(out.maxL, int(p-kLabels));
            }
            char* end=nullptr;
            const long n=std::strtol(h[1].c_str(), &end, 10);
            if (*end || n<=0) return O::Fail("element "+sym+": primitive count '"+h[1]+"'"+where(i));
            std::string scale;
            if (!CNumber(h[2], scale)) return O::Fail("element "+sym+": scale factor '"+h[2]+"'"+where(i));
            if (std::abs(std::strtod(scale.c_str(), nullptr)-1.0)>1e-12)
                return O::Fail("element "+sym+": a shell scale factor of "+h[2]+where(i)+" (the library reads unscaled exponents only)");
            ++i;
            for (long k=0; k<n; ++k, ++i)
            {
                if (i>=lines.size()) return O::Fail("element "+sym+": the "+s.labels+" shell ends before its "+std::to_string(n)+" primitives");
                const auto p=Tokens(lines[i]);
                if (p.size()!=1+s.labels.size())
                    return O::Fail("element "+sym+": a "+s.labels+" primitive needs "+std::to_string(1+s.labels.size())
                                   +" numbers"+where(i)+", found '"+lines[i]+"'");
                std::vector<std::string> row;
                std::string x;
                for (size_t c=0; c<p.size(); ++c)
                {
                    if (!CNumber(p[c], x)) return O::Fail("element "+sym+": '"+p[c]+"' is not a number"+where(i));
                    if (x!=p[c]) hadD=true;
                    (c==0 ? s.exponents.push_back(x) : row.push_back(x));
                }
                s.coeffs.push_back(row);
            }
            shells.push_back(std::move(s));
        }
        if (shells.empty()) return O::Fail("element "+sym+" has no shells");

        // Write the element: one contraction per label (an SP shell is its S and its P on the same exponents -- exact).
        body<<" "<<SYM<<"   0\n";
        for (const Shell& s : shells)
        {
            if (s.labels.size()>1) ++splitSP;
            for (size_t c=0; c<s.labels.size(); ++c)
            {
                body<<" "<<s.labels[c]<<"   "<<s.exponents.size()<<"  1.00\n";
                for (size_t k=0; k<s.exponents.size(); ++k)
                    body<<"   "<<s.exponents[k]<<"   "<<s.coeffs[k][c]<<"\n";
            }
        }
        body<<" ****\n";
    }
    if (out.elements.empty()) return O::Fail("no element block: the text holds no basis");

    if (splitSP) out.notes.push_back(std::to_string(splitSP)+" shared-exponent shell(s) (SP) written as separate S and P contractions");
    if (hadD)    out.notes.push_back("Fortran D exponents written as E");

    std::ostringstream f;
    f<<Preamble(title)<<"\n";
    for (const std::string& p : provenance) f<<Preamble(p)<<"\n";
    for (const std::string& c : comments) if (!c.empty() && c.find_first_not_of('-')!=std::string::npos) f<<c<<"\n";
    f<<"\n!\nBASIS=\""<<Preamble(title)<<"\"\n"<<body.str();
    out.bsd=f.str();
    return O::Ok(std::move(out));
}

Outcome<std::string, std::string> SubsetBasis(const std::string& bsd, const std::vector<std::string>& elements)
{
    using O=Outcome<std::string, std::string>;
    std::set<std::string> want;
    for (std::string e : elements) { for (char& c : e) c=(char)std::toupper((unsigned char)c); want.insert(e); }
    std::istringstream is(bsd);
    std::ostringstream out;
    std::set<std::string> found;
    bool body=false, keep=true, inBlock=false;
    for (std::string l; std::getline(is, l);)
    {
        const auto t=Tokens(l);
        if (!body)                                              // the preamble, through the BASIS= line, verbatim
        {
            out<<l<<"\n";
            if (!t.empty() && t[0].rfind("BASIS=",0)==0) body=true;
            continue;
        }
        if (!inBlock)
        {
            if (t.empty()) continue;
            if (t.size()!=2 || t[1]!="0") break;               // past the last block: a free-text trailer
            std::string sym=t[0];
            for (char& c : sym) c=(char)std::toupper((unsigned char)c);
            keep=want.count(sym)>0;
            if (keep) found.insert(sym);
            inBlock=true;
        }
        if (keep) out<<l<<"\n";
        if (!t.empty() && t[0]=="****") inBlock=false;
    }
    if (!body) return O::Fail("no BASIS= header: not a BasisSetData-layout file");
    for (const std::string& e : want)
        if (!found.count(e)) return O::Fail("the basis has no block for "+e);
    return O::Ok(out.str());
}
}
