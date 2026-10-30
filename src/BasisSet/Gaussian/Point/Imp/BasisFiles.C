// File: BasisSet/Gaussian/Point/Imp/BasisFiles.C
module;
#include <algorithm>
#include <cctype>
#include <cstdio>
#include <fstream>
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
}
