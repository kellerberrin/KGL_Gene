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

private:

  constexpr static const std::string RUF6_FAMILY_{"RUF6"}; // Pf RNA gene family.


};

} // Namespace.

#endif // KGA_ANALYSIS_LIB_UTILITY_H
