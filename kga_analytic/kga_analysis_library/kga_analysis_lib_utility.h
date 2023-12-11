//
// Created by kellerberrin on 27/11/23.
//

#ifndef KGA_ANALYSIS_LIB_UTILITY_H
#define KGA_ANALYSIS_LIB_UTILITY_H

#include "kgl_genome_feature.h"
#include "kgl_genome_genome.h"

namespace kellerberrin::genome::analysis {   //  organization level namespace

// Object cannot be created, just supplies scope and visibility.
// Defines file and time utility functions and misc string/parsing functions.

class KGAUtility {

public:

  KGAUtility() = delete;
  ~KGAUtility() = delete;


  // Return a vector of genes that have a particular text fragment in the description.
  [[nodiscard]] static GeneVector getGeneVector(const std::shared_ptr<const GenomeReference> &genome_ptr,
                                                const std::string &description_text);

  [[nodiscard]] static GeneVector getncRNAGeneVector(const std::shared_ptr<const GenomeReference> &genome_ptr,
                                                     const std::string &desc_uc_text = "",
                                                     size_t max_size = 0);

  [[nodiscard]] static GeneVector getRUF6Genes(const std::shared_ptr<const GenomeReference> &genome_ptr);
  [[nodiscard]] static GeneVector getPFEMP1Genes(const std::shared_ptr<const GenomeReference> &genome_ptr);
  [[nodiscard]] static GeneVector getCircumsporozoite(const std::shared_ptr<const GenomeReference> &genome_ptr);

private:

  constexpr static const std::string RUF6_FAMILY_{"RUF6"}; // Pf ncRNA gene family.
  constexpr static const std::string PFEMP1_FAMILY_{"PFEMP1"}; // Pf var family.
  // PF3D7_0304600, circumsporozoite (CS)
  // Target for RTS,S (37% eff.) and R21 (77% eff.) vaccines.
  // Both of which target the central AA repeat; 'NANP'.
  constexpr static const std::string CIRCUMSPOROZOITE_{"CIRCUMSPORO"};



};

} // Namespace.

#endif // KGA_ANALYSIS_LIB_UTILITY_H
