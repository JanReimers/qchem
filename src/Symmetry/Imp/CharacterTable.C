// File: Symmetry/Imp/CharacterTable.C  The 8 abelian point-group character tables.
module;
#include <string>
#include <vector>
module qchem.Symmetry.CharacterTable;

namespace qchem::Symmetry
{

CharacterTable AbelianCharacterTable(const std::string& s)
{
    // Standard Mulliken labels and +/-1 characters.  Column order matches opTags.
    if (s=="C1")
        return {"C1", {"E"}, {"A"}, {{1}}};

    if (s=="Ci")
        return {"Ci", {"E","i"}, {"Ag","Au"},
                {{1, 1},
                 {1,-1}}};

    if (s=="Cs")
        return {"Cs", {"E","sh"}, {"A'","A''"},
                {{1, 1},
                 {1,-1}}};

    if (s=="C2")
        return {"C2", {"E","C2z"}, {"A","B"},
                {{1, 1},
                 {1,-1}}};

    if (s=="C2h")
        return {"C2h", {"E","C2z","i","sh"}, {"Ag","Bg","Au","Bu"},
                {{1, 1, 1, 1},   // Ag
                 {1,-1, 1,-1},   // Bg
                 {1, 1,-1,-1},   // Au
                 {1,-1,-1, 1}}}; // Bu

    if (s=="C2v")
        return {"C2v", {"E","C2z","sxz","syz"}, {"A1","A2","B1","B2"},
                {{1, 1, 1, 1},   // A1
                 {1, 1,-1,-1},   // A2
                 {1,-1, 1,-1},   // B1
                 {1,-1,-1, 1}}}; // B2

    if (s=="D2")
        return {"D2", {"E","C2z","C2y","C2x"}, {"A","B1","B2","B3"},
                {{1, 1, 1, 1},   // A
                 {1, 1,-1,-1},   // B1
                 {1,-1, 1,-1},   // B2
                 {1,-1,-1, 1}}}; // B3

    if (s=="D2h")
        return {"D2h", {"E","C2z","C2y","C2x","i","sxy","sxz","syz"},
                {"Ag","B1g","B2g","B3g","Au","B1u","B2u","B3u"},
                {{1, 1, 1, 1, 1, 1, 1, 1},   // Ag
                 {1, 1,-1,-1, 1, 1,-1,-1},   // B1g
                 {1,-1, 1,-1, 1,-1, 1,-1},   // B2g
                 {1,-1,-1, 1, 1,-1,-1, 1},   // B3g
                 {1, 1, 1, 1,-1,-1,-1,-1},   // Au
                 {1, 1,-1,-1,-1,-1, 1, 1},   // B1u
                 {1,-1, 1,-1,-1, 1,-1, 1},   // B2u
                 {1,-1,-1, 1,-1, 1, 1,-1}}}; // B3u

    return {"C1", {"E"}, {"A"}, {{1}}};   // fallback: trivial group
}

} //namespace
