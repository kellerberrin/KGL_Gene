//
// Created by kellerberrin on 15/02/18.
//

#ifndef KGL_SEQUENCE_COMPARE_IMPL_H
#define KGL_SEQUENCE_COMPARE_IMPL_H


#include <memory>
#include <string>
#include <map>
#include "kgl_logging.h"
#include "kgl_genome_types.h"
#include "kgl_sequence_base.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


class SequenceComparison {

public:

  explicit SequenceComparison();
  ~SequenceComparison();

  // Comparison

  CompareScore_t MyerHirschbergGlobal(const std::string& sequenceA, const std::string& sequenceB, std::string& compare_str) const;

  CompareScore_t MyerHirschbergLocal(const std::string& sequenceA, const std::string& sequenceB, std::string& compare_str) const;

  CompareScore_t DNALocalAffineGap(const std::string& sequenceA, const std::string& sequenceB, std::string& compare_str) const;

  // Compare sequences in mutation format.

  std::string editItems(const std::string& reference, const std::string& mutant, char delimiter, VariantOutputIndex index_offset) const;

  std::string editDNAItems(std::shared_ptr<const ContigFeatures> contig_ptr,
                           std::shared_ptr<const DNA5SequenceCoding> reference,
                           std::shared_ptr<const DNA5SequenceCoding> mutant,
                           char delimiter,
                           VariantOutputIndex index_offset) const;



private:

  class SequenceManipImpl;       // Forward declaration of the Sequence Manipulation implementation class
  std::unique_ptr<SequenceManipImpl> sequence_manip_impl_ptr_;    // Sequence Manipulation PIMPL


};





}   // namespace genome
}   // namespace kellerberrin









#endif //KGL_KGL_SEQUENCE_COMPARE_IMPL_H
