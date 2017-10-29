//
// Created by kellerberrin on 29/10/17.
//

#include "kgl_mutant_sequence.h"

namespace kgl = kellerberrin::genome;


std::shared_ptr<kgl::DNA5Sequence>
kgl::MutantSequence::mutantSequence(std::shared_ptr<const kgl::DNA5Sequence> sequence_ptr,
                                    const kgl::OffsetVariantMap& variant_map,
                                    const kgl::SortedCDS& sorted_cds) {

  // If no cds then return null sequence.
  if (sorted_cds.empty()) {

    DNA5Sequence::SequenceString null_seq;
    return std::shared_ptr<DNA5Sequence>(std::make_shared<DNA5Sequence>(DNA5Sequence(null_seq)));

  }

// Check bounds.
  if (sorted_cds.rbegin()->second->sequence().end() >= sequence_ptr->length()) {

    ExecEnv::log().error("mutantSequence(), CDS end offset: {} >= target sequence size: {}",
    sorted_cds.rbegin()->second->sequence().end(),
    sequence_ptr->length());
    DNA5Sequence::SequenceString null_seq;
    return std::shared_ptr<DNA5Sequence>(std::make_shared<DNA5Sequence>(DNA5Sequence(null_seq)));

  }

  // For efficiency, pre-allocate storage for the sequence string.
  std::string::size_type calculated_seq_size = 0;
  for (auto cds : sorted_cds) {

    calculated_seq_size += cds.second->sequence().end() - cds.second->sequence().begin();

  }

  DNA5Sequence::SequenceString coding_sequence;
  coding_sequence.reserve(calculated_seq_size + 1); // Just to make sure.

// Get the strand and copy or reverse copy the base complement.
  switch(sorted_cds.begin()->second->sequence().strand()) {

  case StrandSense::UNKNOWN:
  ExecEnv::log().warn("codingSequence(); CDS: {} offset: {} has 'UNKNOWN' ('.') strand assuming 'FORWARD' ('+')",
  sorted_cds.begin()->second->id(),
  sorted_cds.begin()->second->sequence().begin());

  case StrandSense::FORWARD: {

    typename DNA5Sequence::SequenceString::const_iterator begin;
    typename DNA5Sequence::SequenceString::const_iterator end;
    for (auto cds : sorted_cds) {

      begin = sequence_ptr->getSequenceString().begin() + cds.second->sequence().begin();
      end = sequence_ptr->getSequenceString().begin() + cds.second->sequence().end();
      std::copy( begin, end, std::back_inserter(coding_sequence));

    }

  }
  break;

  case StrandSense::REVERSE: {

    // Insert in reverse complement order.
    typename DNA5Sequence::SequenceString::const_reverse_iterator rbegin;
    typename DNA5Sequence::SequenceString::const_reverse_iterator rend;
    auto complement_base =
    [](typename DNA5Sequence::NucleotideType base) { return NucleotideColumn_DNA5::complementNucleotide(base); };
    for (auto rit = sorted_cds.rbegin(); rit != sorted_cds.rend(); ++rit) {

      rbegin = sequence_ptr->getSequenceString().rbegin()
               + (sequence_ptr->getSequenceString().length() - rit->second->sequence().end());
      rend = sequence_ptr->getSequenceString().rbegin()
             + (sequence_ptr->getSequenceString().length() - rit->second->sequence().begin());
      std::transform( rbegin, rend, std::back_inserter(coding_sequence), complement_base);

    }


  }
  break;

  } // switch

  return std::shared_ptr<DNA5Sequence>(std::make_shared<DNA5Sequence>(DNA5Sequence(coding_sequence)));

}
