//
// Created by kellerberrin on 31/10/17.
//


#include "kgl_sequence_base.h"

namespace kgl = kellerberrin::genome;


std::shared_ptr<kgl::DNA5Sequence>
kgl::DNA5Sequence::codingSequence(std::shared_ptr<const DNA5Sequence> base_sequence_ptr,
                                  const SortedCDS& sorted_cds,
                                  ContigOffset_t contig_offset) {

  // If no cds then return null string.
  if (sorted_cds.empty()) {

    SequenceString null_seq;
    return std::shared_ptr<DNA5Sequence>(std::make_shared<DNA5Sequence>(DNA5Sequence(null_seq)));

  }

  // Check bounds.
  if (sorted_cds.rbegin()->second->sequence().end() >= contig_offset + base_sequence_ptr->length()) {

    ExecEnv::log().error("codingSequence(), CDS end offset: {} >= (target sequence size: {} + offset : {})",
                         sorted_cds.rbegin()->second->sequence().end(),
                         base_sequence_ptr->length(),
                         contig_offset);
    SequenceString null_seq;
    return std::shared_ptr<DNA5Sequence>(std::make_shared<DNA5Sequence>(DNA5Sequence(null_seq)));

  }

  // For efficiency, pre-allocate storage for the sequence string.
  std::string::size_type calculated_seq_size = 0;
  for (auto cds : sorted_cds) {

    calculated_seq_size += cds.second->sequence().end() - cds.second->sequence().begin();

  }

  SequenceString coding_sequence;
  coding_sequence.reserve(calculated_seq_size + 1); // Just to make sure.

  // Get the strand and copy or reverse copy the base complement.
  switch(sorted_cds.begin()->second->sequence().strand()) {

    case StrandSense::UNKNOWN:
      ExecEnv::log().warn("codingSequence(); CDS: {} offset: {} has 'UNKNOWN' ('.') strand assuming 'FORWARD' ('+')",
                          sorted_cds.begin()->second->id(),
                          sorted_cds.begin()->second->sequence().begin());

    case StrandSense::FORWARD: {

      typename SequenceString::const_iterator begin;
      typename SequenceString::const_iterator end;
      for (auto cds : sorted_cds) {

        begin = base_sequence_ptr->base_sequence_.begin() + (cds.second->sequence().begin() - contig_offset);
        end = base_sequence_ptr->base_sequence_.begin() + (cds.second->sequence().end() - contig_offset);
        std::copy( begin, end, std::back_inserter(coding_sequence));

      }

    }
      break;

    case StrandSense::REVERSE: {

      // Insert in reverse complement order.
      typename SequenceString::const_reverse_iterator rbegin;
      typename SequenceString::const_reverse_iterator rend;
      auto complement_base =
      [](typename DNA5Sequence::NucleotideType base) { return NucleotideColumn_DNA5::complementNucleotide(base); };
      for (auto rit = sorted_cds.rbegin(); rit != sorted_cds.rend(); ++rit) {

        rbegin = base_sequence_ptr->base_sequence_.rbegin()
                 + (base_sequence_ptr->length() - (rit->second->sequence().end() - contig_offset));
        rend = base_sequence_ptr->base_sequence_.rbegin()
               + (base_sequence_ptr->length() - (rit->second->sequence().begin() - contig_offset));
        std::transform( rbegin, rend, std::back_inserter(coding_sequence), complement_base);

      }


    }
      break;

  } // switch

  // Check the sequence size.
  if (coding_sequence.length() != calculated_seq_size) {

    ExecEnv::log().error("Coding sequence length: {} NOT EQUAL to calculated sequence: {}",
                         coding_sequence.length(),
                         calculated_seq_size);

  }

  return std::shared_ptr<DNA5Sequence>(std::make_shared<DNA5Sequence>(DNA5Sequence(coding_sequence)));

}

// Returns bool false if contig_offset is not within the coding sequences defined by sorted_cds.
// If the contig_offset is in the coding sequence then a valid sequence_offset and the sequence length is returned.
// The offset is adjusted for strand type; the offset arithmetic is reversed for -ve strand sequences.
bool kgl::DNA5Sequence::offsetWithinSequence(const SortedCDS& sorted_cds,
                                             ContigOffset_t contig_offset,
                                             ContigOffset_t& sequence_offset,
                                             ContigSize_t& sequence_length) {

  bool iscoding = false;
  ContigOffset_t coding_offset = 0;
  ContigSize_t coding_size = 0;

  if (sorted_cds.empty()) {

    sequence_offset = 0;
    sequence_length = 0;
    return false;

  }

  // Get the strand.
  StrandSense strand = sorted_cds.begin()->second->sequence().strand();

  switch(strand) {

    case StrandSense::UNKNOWN:  // Complain and assume a forward strand.
    ExecEnv::log().error("CDS feature: {} offset: {} with UNKNOWN strand sense",
                         sorted_cds.begin()->second->id(), sorted_cds.begin()->second->sequence().begin());
    case StrandSense::FORWARD: {

      for (auto cds : sorted_cds) {

        // within the CDS
        if (contig_offset >= cds.second->sequence().begin() and contig_offset <= cds.second->sequence().end()) {

          coding_offset += (contig_offset - cds.second->sequence().begin());
          iscoding = true;

        } else if (contig_offset > cds.second->sequence().end()) {

          coding_offset += (cds.second->sequence().end() - cds.second->sequence().begin());

        }

        coding_size += (cds.second->sequence().end() - cds.second->sequence().begin());

      }

    }
      break;

    case StrandSense::REVERSE: {

      for (auto rit = sorted_cds.rbegin(); rit != sorted_cds.rend(); ++rit) {

        // within the CDS
        if (contig_offset >= rit->second->sequence().begin() and contig_offset <= rit->second->sequence().end()) {

          coding_offset += (rit->second->sequence().end() - contig_offset);
          iscoding = true;

        } else if (contig_offset < rit->second->sequence().begin()) {

          coding_offset += (rit->second->sequence().end() - rit->second->sequence().begin());

        }

        coding_size += (rit->second->sequence().end() - rit->second->sequence().begin());

      }

    }
      break;

  } // switch

  if (iscoding) {

    sequence_offset = coding_offset;
    sequence_length = coding_size;

  } else {

    sequence_offset = 0;
    sequence_length = 0;

  }

  return iscoding;

}


std::shared_ptr<kgl::DNA5Sequence>
kgl::DNA5Sequence::codingSubSequence(std::shared_ptr<const DNA5Sequence> base_sequence_ptr,
                                     const SortedCDS& sorted_cds,
                                     ContigOffset_t& sub_sequence_offset,
                                     ContigSize_t sub_sequence_length,
                                     ContigOffset_t contig_offset) {

  std::shared_ptr<kgl::DNA5Sequence> sequence = codingSequence(base_sequence_ptr, sorted_cds, contig_offset);

  ContigOffset_t sequence_offset;
  ContigSize_t sequence_length;
  if (offsetWithinSequence(sorted_cds, sub_sequence_offset, sequence_offset, sequence_length)) {


  }

  return sequence;

}
