// File: BasisSet/Gaussian/Point/Imp/BasisFiles.C
module;
#include <algorithm>
#include <cctype>
#include <cstdio>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <filesystem>
#ifndef BASISSET_DATA_PATH
#error "BASISSET_DATA_PATH must be defined by CMake"
#endif
module qchem.BasisSet.Gaussian.Point.BasisFiles;

namespace qchem::BasisSet::Gaussian
{
    static const std::filesystem::path theDataDir = BASISSET_DATA_PATH;

    std::string BasisFile(const std::string& filename)
    {
        // A bare name lives in the shipped data directory; anything with a directory part (the deck's \c basis.file
        // override, an absolute or cwd-relative path) is used exactly as given.
        if (filename.find('/')!=std::string::npos) return filename;
        return (theDataDir / filename).string();
    }

    std::vector<std::string> BasisNames()
    {
        std::vector<std::string> names;
        for (const auto& e : std::filesystem::directory_iterator(theDataDir))
            if (e.is_regular_file() && e.path().extension()==".bsd") names.push_back(e.path().stem().string());
        std::sort(names.begin(), names.end());
        return names;
    }

    std::string BasisFileOfName(const std::string& name)
    {
        std::string stem=name;
        for (char& c : stem) c=(char)std::tolower((unsigned char)c);
        if (std::filesystem::is_regular_file(theDataDir/(stem+".bsd"))) return stem+".bsd";
        std::string valid;
        for (const std::string& n : BasisNames()) valid+=(valid.empty()?"":", ")+n;
        throw std::runtime_error("unknown basis \""+name+"\": no "+stem+".bsd in "+theDataDir.string()+"; the names there: "+valid);
    }

    std::string BasisFileChecksum(const std::string& filename)
    {
        const std::string path=BasisFile(filename);
        std::ifstream in(path, std::ios::binary);
        if (!in) throw std::runtime_error("BasisFileChecksum: cannot read basis file "+path);
        unsigned long long h=1469598103934665603ULL; char c;
        while (in.get(c)) { h^=(unsigned char)c; h*=1099511628211ULL; }
        char b[32]; std::snprintf(b, sizeof b, "%016llx", h);
        return b;
    }

    std::string KindName(BasisKind k)
    {
        switch (k)
        {
            case BasisKind::AllElectron: return "all-electron";
            case BasisKind::PPValence:   return "pp-valence";
            case BasisKind::CoulombFit:  return "fit-coulomb";
            case BasisKind::ExchangeFit: return "fit-exchange";
            case BasisKind::Unstated:    break;
        }
        return "";
    }

    BasisFileKind KindOfBasisFile(const std::string& filename)
    {
        const std::string path=BasisFile(filename);
        std::ifstream in(path);
        if (!in) throw std::runtime_error("KindOfBasisFile: cannot read basis file "+path);
        const auto bad=[&](const std::string& why)
        {
            return std::runtime_error("basis file "+path+": "+why+" -- To fix this, please state ONE preamble line "
                "'KIND all-electron', 'KIND pp-valence <pseudopotential set>', 'KIND fit-coulomb' or 'KIND fit-exchange'");
        };
        BasisFileKind out;
        bool seen=false;
        std::string line;
        while (std::getline(in,line))
        {
            std::istringstream is(line);
            std::string tok;
            if (!(is>>tok)) continue;
            if (tok=="!") break;                                   // end of the preamble
            if (tok!="KIND") continue;
            if (seen) throw bad("two KIND lines");
            seen=true;
            std::string word, set, extra;
            is>>word>>set>>extra;
            if (!extra.empty()) throw bad("KIND line '"+line+"' has trailing words");
            for (BasisKind k : {BasisKind::AllElectron, BasisKind::PPValence, BasisKind::CoulombFit, BasisKind::ExchangeFit})
                if (word==KindName(k)) out.kind=k;
            if (out.kind==BasisKind::Unstated) throw bad("unknown KIND '"+word+"'");
            if (out.kind==BasisKind::PPValence && set.empty()) throw bad("KIND pp-valence names no pseudopotential set");
            if (out.kind!=BasisKind::PPValence && !set.empty()) throw bad("KIND "+word+" names a pseudopotential set ("+set+")");
            out.ppSet=set;
        }
        return out;
    }
}
