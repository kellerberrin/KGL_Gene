//
// Created by kellerberrin on 11/12/23.
//

#ifndef KGL_SEQUENCE_MOTIF_H
#define KGL_SEQUENCE_MOTIF_H

#include "kgl_sequence_virtual.h"


namespace kellerberrin::genome{   //  organization level namespace


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// IUPAC nucleotide codes.
//
// A	=> Adenine
// C	=> Cytosine
// G	=> Guanine
// T (or U)	=> Thymine (or Uracil)
// R	=> A or G
// Y	=> C or T
// S	=> G or C
// W	=> A or T
// K	=> G or T
// M	=> A or C
// B	=> C or G or T
// D	=> A or G or T
// H	=> A or C or T
// V	=> A or C or G
// N	=> any base
// . => any base
// - => zero or any base
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// Object cannot be created, just supplies scope and visibility.
class SearchSequence {

public:

  SearchSequence() = delete;
  ~SearchSequence() = delete;

  // Convenience routine to convert IUPAC nucleotide codes into a regex string.
  [[nodiscard]] static std::string IUPACRegex(const std::string_view& IUPAC_search);
  // Search for DNA motifs.
  [[nodiscard]] static std::vector<OpenRightUnsigned> PfPolymerase_III_ABox(const VirtualSequence& sequence) {

    return sequence.regexSearch(IUPACRegex(PF_POL_III_A_BOX_));

  }

  [[nodiscard]] static std::vector<OpenRightUnsigned> PfPolymerase_III_BBox(const VirtualSequence& sequence) {

    return sequence.regexSearch(IUPACRegex(PF_POL_III_B_BOX_));

  }

private:

  constexpr static const std::string PF_POL_III_A_BOX_{"TRGYNNANNNG"}; // Pf Polymerase III 'A' box
  constexpr static const std::string PF_POL_III_B_BOX_{"GWTCRANNC"}; // Pf Polymerase III 'B' box

};


} // Namespace






#endif //KGL_SEQUENCE_MOTIF_H
