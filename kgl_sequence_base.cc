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
                                     ContigOffset_t sub_sequence_offset,  // base count offset; 0 == all
                                     ContigSize_t sub_sequence_length,   // number of bases; 0 == all
                                     ContigOffset_t contig_offset) {

  // If no cds then return null string.
  if (sorted_cds.empty()) {

    SequenceString null_seq;
    return std::shared_ptr<DNA5Sequence>(std::make_shared<DNA5Sequence>(DNA5Sequence(null_seq)));

  }

  // Check bounds.
  if (sorted_cds.rbegin()->second->sequence().end() >= contig_offset + base_sequence_ptr->length()) {

    ExecEnv::log().error("codingSubSequence(), CDS end offset: {} >= (target sequence size: {} + offset : {})",
                         sorted_cds.rbegin()->second->sequence().end(),
                         base_sequence_ptr->length(),
                         contig_offset);
    SequenceString null_seq;
    return std::shared_ptr<DNA5Sequence>(std::make_shared<DNA5Sequence>(DNA5Sequence(null_seq)));

  }

  // Get the size of the coding sequence
  std::string::size_type calculated_seq_size = 0;
  for (auto cds : sorted_cds) {

    calculated_seq_size += cds.second->sequence().end() - cds.second->sequence().begin();

  }

  // If subsequence length is zero and offset is zero then return the whole coding sequence
  if (sub_sequence_offset == 0 and sub_sequence_length == 0) {

    sub_sequence_length = calculated_seq_size;

  }

  // Make sure the requested offset and length are within the coding sequence.
  if ((sub_sequence_offset + sub_sequence_length)  > calculated_seq_size) {

    ExecEnv::log().error("codingSubSequence(), Sub-seq offset: {} + Sub seq length: {} > sequence size: {}",
                         sub_sequence_offset,
                         sub_sequence_length,
                         calculated_seq_size);
    SequenceString null_seq;
    return std::shared_ptr<DNA5Sequence>(std::make_shared<DNA5Sequence>(DNA5Sequence(null_seq)));

  }

  SequenceString coding_sequence;
  coding_sequence.reserve(sub_sequence_length + 1); // Just to make sure.

  // Get the strand and copy or reverse copy the base complement.
  switch(sorted_cds.begin()->second->sequence().strand()) {

    case StrandSense::UNKNOWN:
      ExecEnv::log().warn("codingSubSequence(); CDS: {} offset: {} has 'UNKNOWN' ('.') strand assuming 'FORWARD' ('+')",
                          sorted_cds.begin()->second->id(),
                          sorted_cds.begin()->second->sequence().begin());

    case StrandSense::FORWARD: {

      typename SequenceString::const_iterator begin;
      typename SequenceString::const_iterator end;
      ContigOffset_t begin_offset;
      ContigOffset_t end_offset;
      ContigOffset_t relative_offset = 0;
      // Convert to an absolute sequence based offset
      for (auto cds : sorted_cds) {

        ContigSize_t cds_size = cds.second->sequence().end() - cds.second->sequence().begin();

        if (sub_sequence_offset < relative_offset + cds_size) {

          if (sub_sequence_offset <= relative_offset) {

            begin_offset = cds.second->sequence().begin();

          } else {

            begin_offset = cds.second->sequence().begin() + (sub_sequence_offset - relative_offset);

          }

          if (sub_sequence_offset + sub_sequence_length > relative_offset + cds_size) {

            end_offset = cds.second->sequence().end();

          } else {

            end_offset = cds.second->sequence().end()
                         - ((relative_offset + cds_size) - (sub_sequence_offset + sub_sequence_length));

          }

          begin = base_sequence_ptr->base_sequence_.begin() + (begin_offset - contig_offset);
          end = base_sequence_ptr->base_sequence_.begin() + (end_offset - contig_offset);
          std::copy( begin, end, std::back_inserter(coding_sequence));

        } // if sub_sequence_offset < relative_offset + cds_size

        relative_offset += cds_size;
        // Check if we need to process more CDS.
        if (sub_sequence_offset + sub_sequence_length < relative_offset) {

          break;

        } // terminate if complete.

      } // for cds

    } // case
      break;

    case StrandSense::REVERSE: {

      // Insert in reverse complement order.
      typename SequenceString::const_reverse_iterator rbegin;
      typename SequenceString::const_reverse_iterator rend;
      ContigOffset_t begin_offset;
      ContigOffset_t end_offset;
      ContigOffset_t relative_offset = 0;
      auto complement_base =
      [](typename DNA5Sequence::NucleotideType base) { return NucleotideColumn_DNA5::complementNucleotide(base); };
      for (auto rit = sorted_cds.rbegin(); rit != sorted_cds.rend(); ++rit) {

        ContigSize_t cds_size = rit->second->sequence().end() - rit->second->sequence().begin();

        if (sub_sequence_offset < relative_offset + cds_size) {

          if (sub_sequence_offset <= relative_offset) {

            begin_offset = rit->second->sequence().end();

          } else {

            begin_offset = rit->second->sequence().end() - (sub_sequence_offset - relative_offset);

          }

          if (sub_sequence_offset + sub_sequence_length > relative_offset + cds_size) {

            end_offset = rit->second->sequence().begin();

          } else {

            end_offset = rit->second->sequence().begin()
                         + ((relative_offset + cds_size) - (sub_sequence_offset + sub_sequence_length));

          }

          rbegin = base_sequence_ptr->base_sequence_.rbegin()
                   + (base_sequence_ptr->length() - (begin_offset - contig_offset));
          rend = base_sequence_ptr->base_sequence_.rbegin()
                 + (base_sequence_ptr->length() - (end_offset - contig_offset));
          std::transform( rbegin, rend, std::back_inserter(coding_sequence), complement_base);

        } // if sub_sequence_offset < relative_offset + cds_size

        relative_offset += cds_size;
        // Check if we need to process more CDS.
        if (sub_sequence_offset + sub_sequence_length < relative_offset) {

          break;

        } // terminate if complete

      } // for cds

    } // case
      break;

  } // switch

  // Check the sub sequence size.
  if (coding_sequence.length() != sub_sequence_length) {

    ExecEnv::log().error("Coding SubSequence length: {} NOT EQUAL to specified subsequence: {}",
                         coding_sequence.length(),
                         calculated_seq_size);

  }

  return std::shared_ptr<DNA5Sequence>(std::make_shared<DNA5Sequence>(DNA5Sequence(coding_sequence)));

}


bool kgl::DNA5Sequence::modifyBase(NucleotideType Nucleotide, ContigOffset_t sequence_offset) {

  if (sequence_offset >= base_sequence_.length()) {

    ExecEnv::log().error("modifyBase(), sequence offset: {} exceeds sequence size: {}",
                         sequence_offset, base_sequence_.length());
    return false;
  }

  base_sequence_.at(sequence_offset) = Nucleotide;
  return true;

}
