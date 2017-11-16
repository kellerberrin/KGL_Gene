//
// Created by kellerberrin on 16/11/17.
//



#include "kgl_sequence_base_alpha.h"
#include "kgl_sequence_codon_alpha.h"

namespace kgl = kellerberrin::genome;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The base DNA5 sequence class.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool kgl::DNA5Sequence::modifyBase(DNA5::Alphabet nucleotide, ContigOffset_t sequence_offset) {

  if (sequence_offset >= base_sequence_.length()) {

    ExecEnv::log().error("modifyBase(), sequence offset: {} exceeds sequence size: {}",
                         sequence_offset, base_sequence_.length());
    return false;
  }

  base_sequence_.modifyNucleotide(sequence_offset, nucleotide);

  return true;

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// A linear and contiguous DNA5 sequence that cannot be used to directly generate an Amino Acid sequence
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////



std::shared_ptr<kgl::DNA5SequenceCoding>
kgl::DNA5SequenceLinear::codingSubSequence(std::shared_ptr<const DNA5SequenceLinear> base_sequence_ptr,
                                           std::shared_ptr<const CodingSequence> coding_seq_ptr,
                                           ContigOffset_t sub_sequence_offset,  // base count offset; 0 == all
                                           ContigSize_t sub_sequence_length,   // number of bases; 0 == all
                                           ContigOffset_t contig_offset) {


  const SortedCDS& sorted_cds = coding_seq_ptr->getSortedCDS();
  // If no cds then return null string.
  if (sorted_cds.empty()) {

    StringDNA5 null_seq;
    return std::shared_ptr<DNA5SequenceCoding>(std::make_shared<DNA5SequenceCoding>(DNA5SequenceCoding(null_seq)));

  }

  // Check bounds.
  if (sorted_cds.rbegin()->second->sequence().end() >= contig_offset + base_sequence_ptr->length()) {

    ExecEnv::log().error("codingSubSequence(), CDS end offset: {} >= (target sequence size: {} + offset : {})",
                         sorted_cds.rbegin()->second->sequence().end(),
                         base_sequence_ptr->length(),
                         contig_offset);
    StringDNA5 null_seq;
    return std::shared_ptr<DNA5SequenceCoding>(std::make_shared<DNA5SequenceCoding>(DNA5SequenceCoding(null_seq)));

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
    StringDNA5 null_seq;
    return std::shared_ptr<DNA5SequenceCoding>(std::make_shared<DNA5SequenceCoding>(DNA5SequenceCoding(null_seq)));

  }

  StringDNA5 coding_sequence;
  coding_sequence.reserve(sub_sequence_length + 1); // Just to make sure.

  // Get the strand and copy or reverse copy the base complement.
  switch(coding_seq_ptr->getGene()->sequence().strand()) {

    case StrandSense::UNKNOWN:
      ExecEnv::log().warn("codingSubSequence(); Gene: {} offset: {} has 'UNKNOWN' ('.') strand assuming 'FORWARD' ('+')",
                          coding_seq_ptr->getGene()->id(),
                          coding_seq_ptr->getGene()->sequence().begin());

    case StrandSense::FORWARD: {

      StringDNA5::const_iterator begin;
      StringDNA5::const_iterator end;
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
      StringDNA5::const_reverse_iterator rbegin;
      StringDNA5::const_reverse_iterator rend;
      ContigOffset_t begin_offset;
      ContigOffset_t end_offset;
      ContigOffset_t relative_offset = 0;

      auto complement_base = [](DNA5::Alphabet base) { return DNA5::complementNucleotide(base); };

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

  return std::shared_ptr<DNA5SequenceCoding>(std::make_shared<DNA5SequenceCoding>(DNA5SequenceCoding(coding_sequence)));

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// A linear and contiguous DNA5 sequence used in a contig (chromosome). This object exists for semantic reasons.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// Returns the codon offset of offset within a coding, returns false if not within the coding sequence.
bool kgl::DNA5SequenceContig::codonOffset(std::shared_ptr<const CodingSequence> coding_seq_ptr,
                                          ContigOffset_t contig_offset,
                                          ContigOffset_t& codon_offset,
                                          ContigSize_t& base_in_codon) const {

  ContigOffset_t sequence_offset;
  ContigSize_t sequence_length;
  if (offsetWithinSequence(coding_seq_ptr, contig_offset, sequence_offset, sequence_length)) {

    codon_offset = static_cast<ContigOffset_t>(sequence_offset / Codon::CODON_SIZE);
    base_in_codon = static_cast <ContigOffset_t>(sequence_offset % Codon::CODON_SIZE);
    return true;

  } else {

    codon_offset = 0;
    base_in_codon = 0;
    return false;

  }

}

// Returns bool false if contig_offset is not within the coding sequence defined by the coding_seq_ptr.
// If the contig_offset is in the coding sequence then a valid sequence_offset and the sequence length is returned.
// The offset is adjusted for strand type; the offset arithmetic is reversed for -ve strand sequences.
bool kgl::DNA5SequenceContig::offsetWithinSequence(std::shared_ptr<const CodingSequence> coding_seq_ptr,
                                                   ContigOffset_t contig_offset,
                                                   ContigOffset_t& sequence_offset,
                                                   ContigSize_t& sequence_length) const {

  bool iscoding = false;
  ContigOffset_t coding_offset = 0;
  ContigSize_t coding_size = 0;
  const SortedCDS sorted_cds = coding_seq_ptr->getSortedCDS();

  if (sorted_cds.empty()) {

    sequence_offset = 0;
    sequence_length = 0;
    return false;

  }


  // Get the strand.
  StrandSense strand = coding_seq_ptr->getGene()->sequence().strand();

  switch(strand) {

    case StrandSense::UNKNOWN:  // Complain and assume a forward strand.
      ExecEnv::log().error("Gene feature: {} offset: {} with UNKNOWN strand sense",
                           coding_seq_ptr->getGene()->id(), coding_seq_ptr->getGene()->sequence().begin());
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

      // Careful with this logic as the CDS offsets are [begin, end). Therefore we must adjust the offset by -1.
    case StrandSense::REVERSE: {

      for (auto rit = sorted_cds.rbegin(); rit != sorted_cds.rend(); ++rit) {

        // within the CDS
        if (contig_offset >= rit->second->sequence().begin() and contig_offset < rit->second->sequence().end()) {

          coding_offset += (rit->second->sequence().end() - contig_offset) - 1;  // Careful here.
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
