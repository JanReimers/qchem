// File: BasisSet/Gaussian/Point/Readers/Gaussian94Import.C  A Gaussian94 basis from OUTSIDE (the Basis Set Exchange's
// "gaussian94" format, a journal's supplement) translated into the BasisSetData layout the Gaussian94Reader reads.
//
// OpenWork §2 "GUI: browse and download a basis from the Basis Set Exchange" (user 2026-11-05; lane G 2026-11-10).  The
// Exchange serves the same layout our .bsd files use, but not exactly: Fortran D exponents (1.30D+01), element symbols in
// mixed case (Fe), a header of "!" comment lines where the reader expects ONE "!" and a BASIS= title, and Pople SP shells
// whose S and P coefficients differ (the reader builds a multi-L shell only when they agree).  The importer is a TRANSLATOR:
// it rewrites those four, and refuses by name what the library cannot represent -- an ECP block (there is no ECP route), a
// shell scale factor other than 1, an angular label beyond I -- instead of dropping it.  An SP shell is split into its S and
// P contractions on the same exponents, which is the same function space exactly.
//
// The input is UNTRUSTED TEXT (a network download): it is parsed as data, never executed, and nothing in it chooses a path.
module;
#include <string>
#include <vector>
export module qchem.BasisSet.Gaussian.Point.Readers.Gaussian94Import;
export import qchem.Outcome;

export namespace qchem::BasisSet::Gaussian
{
//! One imported basis: the file text and what it holds.
struct ImportedBasis
{
    std::string              bsd;        //!< the .bsd text: title and provenance preamble, "!", BASIS="<title>", element blocks
    std::vector<std::string> elements;   //!< element symbols in file order, as the reader spells them (upper case)
    int                      maxL = 0;   //!< the highest angular momentum in the set (0 = s)
    std::vector<std::string> notes;      //!< what the translation rewrote, one line each (SP shells split, D exponents)
};

//! \brief TRANSLATE Gaussian94 text into the BasisSetData layout.  \a title becomes the BASIS= line and the first line;
//! \a provenance lines (source, version, citations) go into the free-text preamble the reader skips, as do the input's own
//! "!" comment lines.  FAILS (the reason, naming the element and shell) on an ECP block, a scale factor other than 1, an
//! unknown angular label, a malformed number or count, a repeated element, or no element at all.
Outcome<ImportedBasis, std::string> ImportGaussian94(const std::string& text, const std::string& title,
                                                     const std::vector<std::string>& provenance = {});

//! \brief The part of a BasisSetData-layout file a run USES: the whole preamble (title, attribution, PP-PROVENANCE, the "!" and
//! BASIS= lines) verbatim and the blocks of \a elements (symbols, any case), in file order -- what a deck inlines (user ruling
//! 2026-11-10: a deck carries a non-shipped basis inline, for the elements it uses, with its attribution).  FAILS naming an element
//! the file has no block for.
Outcome<std::string, std::string> SubsetBasis(const std::string& bsd, const std::vector<std::string>& elements);
}
