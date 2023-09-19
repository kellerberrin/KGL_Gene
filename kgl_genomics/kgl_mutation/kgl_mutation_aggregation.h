//
// Created by kellerberrin on 16/09/23.
//

#ifndef KGL_MUTATION_AGGREGATION_H
#define KGL_MUTATION_AGGREGATION_H

#include "kgl_mutation_sequence.h"
#include "kgl_genome_interval.h"


namespace kellerberrin::genome {   //  organization::project level namespace


class SequenceAggregation {


public:

  SequenceAggregation() = default;
  ~SequenceAggregation() = default;

  // Given a vector of intervals, sorted them by lower() and then retrieves the modified sequence,
  // and concatenates these in sorted sequence order.
  static std::optional<DNA5SequenceLinear> concatModifiedSequences( const AdjustedSequence& adjusted_sequence,
                                                                    const std::vector<OpenRightUnsigned>& interval_vector);

  // Given a vector of intervals, sorted them by lower() and then retrieves the original sequence,
  // and concatenates these in sorted sequence order.
  static std::optional<DNA5SequenceLinear> concatOriginalSequences( const AdjustedSequence& adjusted_sequence,
                                                                    const std::vector<OpenRightUnsigned>& interval_vector);

  static std::optional<DNA5SequenceLinear> getModifiedGene( const AdjustedSequence& adjusted_sequence,
                                                            const GeneIntervalStructure& gene_interval,
                                                            const FeatureIdent_t& transcript_id);

  static std::optional<DNA5SequenceLinear> getOriginalGene( const AdjustedSequence& adjusted_sequence,
                                                            const GeneIntervalStructure& gene_interval,
                                                            const FeatureIdent_t& transcript_id);

private:


};



} // namespace

#endif //KGL_MUTATION_AGGREGATION_H
