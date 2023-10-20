//
// Created by kellerberrin on 20/10/23.
//

#ifndef KGL_GENOME_VERIFY_H
#define KGL_GENOME_VERIFY_H

#include "kgl_genome_prelim.h"

#include <vector>

namespace kellerberrin::genome {   //  organization level namespace



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Simple object to hold coding sequence validity statistics.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct SequenceValidityStatistics {

public:

  SequenceValidityStatistics() =default;
  ~SequenceValidityStatistics() =default;

  [[nodiscard]] size_t ncRNA() const { return ncRNA_; }
  [[nodiscard]] size_t validProtein() const { return valid_protein_; }
  [[nodiscard]] size_t empty() const { return empty_; }
  [[nodiscard]] size_t notMod3() const { return not_mod3_; }
  [[nodiscard]] size_t noStartCodon() const { return no_start_codon_; }
  [[nodiscard]] size_t nonsenseMutation() const { return nonsense_mutation_; }
  [[nodiscard]] const std::vector<size_t>& nonsenseMutationSize() const { return nonsense_mutation_size_; }
  [[nodiscard]] size_t noStopCodon() const { return no_stop_codon_; }
  [[nodiscard]] size_t invalidProtein() const { return noStartCodon() + empty() + notMod3() + nonsenseMutation() +
                                                       noStopCodon(); }
  [[nodiscard]] size_t totalProtein() const { return validProtein() + invalidProtein(); }
  [[nodiscard]] size_t totalSequence() const { return totalProtein() + ncRNA(); }

  void updateValidity(CodingSequenceValidity validity, size_t amino_size = 0);
  void updateTranscriptArray(const std::shared_ptr<const TranscriptionSequenceArray>& transcript_array_ptr);
  void updateTranscript(const std::shared_ptr<const TranscriptionSequence>& transcript_ptr);

private:


  size_t ncRNA_{0};
  size_t valid_protein_{0};
  size_t empty_{0};
  size_t not_mod3_{0};
  size_t no_start_codon_{0};
  size_t nonsense_mutation_{0};
  std::vector<size_t> nonsense_mutation_size_;
  size_t no_stop_codon_{0};


};



} // Namespace

#endif //KGL_GENOME_VERIFY_H
