//
// Created by kellerberrin on 7/09/23.
//

#ifndef KGL_MUTATION_SEQUENCE_H
#define KGL_MUTATION_SEQUENCE_H

#include "kgl_mutation_variant_map.h"
#include "kgl_genome_contig.h"


namespace kellerberrin::genome {   //  organization::project level namespace


class AdjustedSequence {

public:

  AdjustedSequence( const std::shared_ptr<const ContigReference>& contig_ref_ptr,
                    OpenRightInterval contig_interval,
                    IntervalModifyMap interval_modify_map)
                    : contig_interval_(contig_interval), interval_modify_map_(std::move(interval_modify_map)) {

    if (not interval_modify_map_.empty()) {

      auto [offset, first_sequence_modify] = *interval_modify_map_.begin();

      if (contig_interval != first_sequence_modify.priorInterval()) {

        ExecEnv::log().error("AdjustedSequence::AdjustedSequence; contig interval: {} does not match initial map interval: {}",
                             contig_interval.toString(), first_sequence_modify.priorInterval().toString());
        // Clear the map for good measure.
        interval_modify_map_.clear();

      }

    }

    modified_sequence_ = contig_ref_ptr->getSubSequence(contig_interval);

  }
  ~AdjustedSequence() = default;

  [[nodiscard]] OpenRightInterval contigInterval() const { return contig_interval_; }
  [[nodiscard]] const IntervalModifyMap& intervalMap() const { return interval_modify_map_; }
  [[nodiscard]] const DNA5SequenceLinear& modifiedSequence() const { return modified_sequence_; }

  [[nodiscard]] bool updateSequence();

private:


  OpenRightInterval contig_interval_;
  IntervalModifyMap interval_modify_map_;
  DNA5SequenceLinear modified_sequence_;

};



}  // Namespace



#endif //KGL_MUTATION_SEQUENCE_H
