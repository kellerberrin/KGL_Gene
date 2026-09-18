//
// kgl_sequence_motif.h — IUPAC nucleotide code to regex conversion and Pf motif search.
//
// SearchSequence is retained as a namespace so the reference call spelling is unchanged.
//

#ifndef KGL_SEQUENCE_MOTIF_H
#define KGL_SEQUENCE_MOTIF_H


#include <regex>
#include <string>
#include <string_view>
#include <vector>

#include "kgl_genome_types.h"
#include "kgl_sequence_view.h"
#include "kel_interval_unsigned.h"


namespace kellerberrin::genome {   //  organization::project level namespace


// A namespace (not a deleted-constructor class): the call sites
// SearchSequence::PfPolymerase_III_ABox(...) compile unchanged.
namespace SearchSequence {

  inline constexpr std::string_view PF_POL_III_A_BOX_{"TRGYNNANNNG"};   // Pf Polymerase III 'A' box
  inline constexpr std::string_view PF_POL_III_B_BOX_{"GWTCRANNC"};    // Pf Polymerase III 'B' box

  /// Convenience routine to convert IUPAC nucleotide codes into a regex string.
  [[nodiscard]] std::string IUPACRegex(std::string_view IUPAC_search);

  /// Search for DNA motifs.
  [[nodiscard]] inline std::vector<OpenRightUnsigned> PfPolymerase_III_ABox(SequenceLike auto const& sequence) {

    static const std::regex compiled{IUPACRegex(PF_POL_III_A_BOX_)};
    return regexSearch(sequence, compiled);

  }

  /// Search for DNA motifs.
  [[nodiscard]] inline std::vector<OpenRightUnsigned> PfPolymerase_III_BBox(SequenceLike auto const& sequence) {

    static const std::regex compiled{IUPACRegex(PF_POL_III_B_BOX_)};
    return regexSearch(sequence, compiled);

  }

}   // namespace SearchSequence


}   // end namespace


#endif //KGL_SEQUENCE_MOTIF_H
