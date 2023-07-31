//
// Created by kellerberrin on 3/01/18.
//

#include "kgl_sequence_offset.h"

namespace kgl = kellerberrin::genome;


// Returns a defined subsequence (generally a single/group of codons) of the coding sequence
// Setting sub_sequence_offset and sub_sequence_length to zero copies the entire sequence defined by the TranscriptionFeatureMap.
kgl::DNA5SequenceCoding kgl::SequenceOffset::refCodingSubSequence( const std::shared_ptr<const TranscriptionSequence>& coding_seq_ptr,
                                                                   const DNA5SequenceLinear& sequence_ptr,
                                                                   ContigOffset_t sub_sequence_offset,
                                                                   ContigSize_t sub_sequence_length,
                                                                   ContigOffset_t contig_offset) {

  StrandSense strand;
  ExonOffsetMap exon_offset_map;
  if (not exonOffsetAdapter(coding_seq_ptr, strand, exon_offset_map)) {

    ExecEnv::log().error("SequenceOffset::refCodingSubSequence; DNA Sequence length: {} cannot generate EXON map for subsequence offset: {}, length: {}",
                         sequence_ptr.length(), sub_sequence_offset, sub_sequence_length);
    exon_offset_map.clear();

  }

  return codingSubSequence(sequence_ptr, exon_offset_map, coding_seq_ptr->strand(), sub_sequence_offset, sub_sequence_length, contig_offset);

}

// Convenience routine that returns an array of sequences (strand adjusted) from an unmutated (reference) sequence.
// Returned sequences are in transcription (strand) order with array[0] being the first exon.
 std::vector<kgl::DNA5SequenceCoding> kgl::SequenceOffset::refExonArraySequence( const std::shared_ptr<const TranscriptionSequence>& coding_seq_ptr,
                                                                                 const DNA5SequenceLinear& sequence_ptr,
                                                                                 ContigOffset_t sub_sequence_offset,
                                                                                 ContigSize_t sub_sequence_length,
                                                                                 ContigOffset_t contig_offset) {

  std::vector<kgl::DNA5SequenceCoding> exon_vector;
  StrandSense strand;
  ExonOffsetMap exon_offset_map;

  if (not exonOffsetAdapter(coding_seq_ptr, strand, exon_offset_map)) {

    ExecEnv::log().error("SequenceOffset::refExonArraySequence; DNA Sequence length: {} cannot generate EXON map for subsequence offset: {}, length: {}",
                         sequence_ptr.length(), sub_sequence_offset, sub_sequence_length);
    exon_offset_map.clear();

  }

  for (auto exon_offset : exon_offset_map) {

    ExonOffsetMap offset_map = { exon_offset };
    DNA5SequenceCoding exon_sequence = codingSubSequence( sequence_ptr,
                                                          offset_map,
                                                          coding_seq_ptr->strand(),
                                                          sub_sequence_offset,
                                                          sub_sequence_length,
                                                          contig_offset);
    exon_vector.push_back(std::move(exon_sequence));

  }

  return exon_vector;

}



// Returns a defined subsequence of all the introns of the coding sequence
// Setting sub_sequence_offset and sub_sequence_length to zero copies the entire intron sequence defined by the TranscriptionFeatureMap.
kgl::DNA5SequenceCoding kgl::SequenceOffset::refIntronSubSequence( const std::shared_ptr<const TranscriptionSequence>& coding_seq_ptr,
                                                                   const DNA5SequenceLinear& sequence_ptr,
                                                                   ContigOffset_t sub_sequence_offset,
                                                                   ContigSize_t sub_sequence_length,
                                                                   ContigOffset_t contig_offset) {

  StrandSense strand;
  IntronOffsetMap intron_offset_map;
  if (not intronOffsetAdapter(coding_seq_ptr, strand, intron_offset_map)) {

    ExecEnv::log().error("refIntronSubSequence(), DNA Sequence length: {} cannot generate INTRON map for subsequence offset: {}, length: {}",
                         sequence_ptr.length(), sub_sequence_offset, sub_sequence_length);
    intron_offset_map.clear();

  }

  return codingSubSequence(sequence_ptr, intron_offset_map, coding_seq_ptr->strand(), sub_sequence_offset, sub_sequence_length, contig_offset);

}


// Convenience routine that returns an array of intron (strand adjusted) from an unmutated (reference) sequence.
// Returned sequences are in transcription (strand) order with array[0] being the first intron.
std::vector<kgl::DNA5SequenceCoding> kgl::SequenceOffset::refIntronArraySequence( const std::shared_ptr<const TranscriptionSequence>& coding_seq_ptr,
                                                                                  const DNA5SequenceLinear& sequence_ptr,
                                                                                  ContigOffset_t sub_sequence_offset,
                                                                                  ContigSize_t sub_sequence_length,
                                                                                  ContigOffset_t contig_offset) {

  std::vector<kgl::DNA5SequenceCoding> intron_vector;
  StrandSense strand;
  IntronOffsetMap intron_offset_map;
  if (not intronOffsetAdapter(coding_seq_ptr, strand, intron_offset_map)) {

    ExecEnv::log().error("SequenceOffset::refIntronArraySequence; DNA Sequence length: {} cannot generate INTRON map for subsequence offset: {}, length: {}",
                         sequence_ptr.length(), sub_sequence_offset, sub_sequence_length);
    intron_offset_map.clear();

  }

  for (auto intron_offset : intron_offset_map) {

    ExonOffsetMap offset_map = { intron_offset };
    DNA5SequenceCoding exon_sequence = codingSubSequence( sequence_ptr,
                                                          offset_map,
                                                          coding_seq_ptr->strand(),
                                                          sub_sequence_offset,
                                                          sub_sequence_length,
                                                          contig_offset);
    intron_vector.push_back(std::move(exon_sequence));

  }

  return intron_vector;

}


kgl::DNA5SequenceCoding kgl::SequenceOffset::mutantCodingSubSequence( const std::shared_ptr<const TranscriptionSequence>& coding_seq_ptr,
                                                                      const DNA5SequenceLinear& sequence_ptr,
                                                                      const VariantMutationOffset& indel_adjust,
                                                                      ContigOffset_t sub_sequence_offset,
                                                                      ContigSize_t sub_sequence_length,
                                                                      ContigOffset_t contig_offset) {

  StrandSense strand;
  ExonOffsetMap exon_offset_map;
  if (not exonMutantOffset(coding_seq_ptr, indel_adjust, strand, exon_offset_map)) {

    ExecEnv::log().error("SequenceOffset::mutantCodingSubSequence; DNA Sequence length: {} cannot generate EXON map for INDEL adjusted subsequence offset: {}, length: {}",
                         sequence_ptr.length(), sub_sequence_offset, sub_sequence_length);
    exon_offset_map.clear();

  }

  return codingSubSequence(sequence_ptr, exon_offset_map, coding_seq_ptr->strand(), sub_sequence_offset, sub_sequence_length, contig_offset);

}



// A sorted map of coding sequence offsets.
// The first (sort key) offset is the beginning of the CDS and the second (value) is the CDS end offset.
// The offsets use the half interval idiom [start, end).
// The offsets are initially contig offsets. However these can be modified to account for indel offset modification.
// The map is generated by an adapter function exonOffsetAdapter() that uses a coding sequence as input.
bool kgl::SequenceOffset::exonOffsetAdapter( const std::shared_ptr<const TranscriptionSequence>& coding_seq_ptr,
                                             StrandSense& strand,
                                             ExonOffsetMap& exon_offset_map) {

  bool return_result = true;

  strand = coding_seq_ptr->strand();


  // All CDS offsets are relative to the underlying contig and are half intervals [begin, end).
  for (auto const& [cds_offset, cds_feature_ptr] : coding_seq_ptr->getFeatureMap()) {

    std::pair<ContigOffset_t, ContigOffset_t> exon_offsets(cds_feature_ptr->sequence().begin(), cds_feature_ptr->sequence().end());

    auto [iter, result] = exon_offset_map.insert(exon_offsets);

    if (not result) {

      ExecEnv::log().error("SequenceOffset::exonOffsetAdapter; Duplicate exon offset: {} for coding sequence id:",
                           exon_offsets.first, coding_seq_ptr->getParent()->id());
      return_result = false;

    }

  }

  return return_result;

}


// A sorted map of intron coding sequence offsets.
// The first (sort key) offset is the beginning of the intron and the second (value) is the intron end offset.
// The offsets use the half interval idiom [start, end).
// The offsets are initially contig offsets. However these can be modified to account for indel offset modification.
// The map is generated by an adapter function intronOffsetAdapter() that uses a coding sequence as input.
// If the coding sequence only defines 1 exon then an empty map is returned.
// In general the number of intron offset records is the number of defined sequences - 1
bool kgl::SequenceOffset::intronOffsetAdapter( const std::shared_ptr<const TranscriptionSequence>& coding_seq_ptr,
                                               StrandSense& strand,
                                               IntronOffsetMap& intron_offset_map) {

  bool return_result = true;

  strand = coding_seq_ptr->strand();

  intron_offset_map.clear();

  // If exon size <= 1 then no introns.
  if (coding_seq_ptr->getFeatureMap().size() <= 1) {

    return return_result;

  }

  auto CDS_iter = coding_seq_ptr->getFeatureMap().begin();

  // All CDS offsets are relative to the underlying contig and are half intervals [begin, end).
  while (CDS_iter != coding_seq_ptr->getFeatureMap().end()) {

    ContigOffset_t begin_intron = CDS_iter->second->sequence().end();

    ++CDS_iter;

    if (CDS_iter == coding_seq_ptr->getFeatureMap().end()) {

      continue;

    }

    ContigOffset_t end_intron = CDS_iter->second->sequence().begin();

    std::pair<ContigOffset_t, ContigOffset_t> intron_offsets(begin_intron, end_intron);

    auto result = intron_offset_map.insert(intron_offsets);

    if (not result.second) {

      ExecEnv::log().error("SequenceOffset::intronOffsetAdapter; Duplicate intron offset: {} for coding sequence id:",
                           intron_offsets.first, coding_seq_ptr->getParent()->id());
      return_result = false;

    }

    ExecEnv::log().info("SequenceOffset::intronOffsetAdapter; Intron: {} [{}, {}) ({})",
                        coding_seq_ptr->getParent()->id(), intron_offsets.first, intron_offsets.second,
                        intron_offsets.second - intron_offsets.first);

  }

  ExecEnv::log().info("SequenceOffset::intronOffsetAdapter; Introns: {}", intron_offset_map.size());

  return return_result;

}



// Adjusted for indel mutations.
bool kgl::SequenceOffset::exonMutantOffset(const std::shared_ptr<const TranscriptionSequence>& coding_seq_ptr,
                                           const VariantMutationOffset& indel_adjust,
                                           StrandSense& strand,
                                           ExonOffsetMap& exon_offset_map) {

  bool return_result = true;

  ExonOffsetMap ref_exon_offset_map;
  return_result = exonOffsetAdapter(coding_seq_ptr, strand, ref_exon_offset_map);

  for (auto exon : ref_exon_offset_map) {

    ContigOffset_t begin_offset = exon.first + indel_adjust.adjustIndelOffsets(exon.first);
    ContigOffset_t end_offset = exon.second + indel_adjust.adjustIndelOffsets(exon.second);

    auto result = exon_offset_map.insert(std::pair<ContigOffset_t , ContigOffset_t >(begin_offset, end_offset));

    if (not result.second) {

      ExecEnv::log().error("exonMutantOffset(), Duplicate exon offset: {} for coding sequence id:",
                           begin_offset, coding_seq_ptr->getParent()->id());
      return_result = false;

    }

  }

  return return_result;

}


kgl::DNA5SequenceCoding kgl::SequenceOffset::codingSubSequence(const DNA5SequenceLinear& base_sequence,
                                                               const ExonOffsetMap& exon_offset_map,
                                                               StrandSense strand,
                                                               ContigOffset_t sub_sequence_offset,  // base count offset; 0 == all
                                                               ContigSize_t sub_sequence_length,   // number of bases; 0 == all
                                                               ContigOffset_t contig_offset) {


  // If no cds then return null string.
  if (exon_offset_map.empty()) {

    return DNA5SequenceCoding();

  }

  // Check bounds.
  if (exon_offset_map.rbegin()->second > contig_offset + base_sequence.length()) {

    ExecEnv::log().error("SequenceOffset::codingSubSequence(), CDS end offset: {} > (target sequence size: {} + offset : {})",
                         exon_offset_map.rbegin()->second,
                         base_sequence.length(),
                         contig_offset);

    return DNA5SequenceCoding();

  }

  // Get the size of the coding sequence
  std::string::size_type calculated_seq_size = 0;
  for (auto exon : exon_offset_map) {

    calculated_seq_size += exon.second - exon.first;

  }

  // If subsequence length is zero and offset is zero then return the whole coding sequence
  if (sub_sequence_offset == 0 and sub_sequence_length == 0) {

    sub_sequence_length = calculated_seq_size;

  }

  // Make sure the requested offset and length are within the coding sequence.
  if ((sub_sequence_offset + sub_sequence_length)  > calculated_seq_size) {

    ExecEnv::log().error("SequenceOffset::codingSubSequence(); Sub-seq offset: {} + Sub seq length: {} > sequence size: {}",
                         sub_sequence_offset,
                         sub_sequence_length,
                         calculated_seq_size);

    return DNA5SequenceCoding();

  }

  StringCodingDNA5 coding_sequence;
  coding_sequence.reserve(sub_sequence_length + 1); // Just to make sure.


  // Get the strand and copy or reverse copy the base complement.
  switch(strand) {

    case StrandSense::FORWARD: {

      StringDNA5::const_iterator begin;
      StringDNA5::const_iterator end;
      ContigOffset_t begin_offset;
      ContigOffset_t end_offset;
      ContigOffset_t relative_offset = 0;

      auto convert_base = [](DNA5::Alphabet base) { return DNA5::convertToCodingDNA5(base); };

      // Convert to an absolute sequence based offset
      for (auto exon : exon_offset_map) {

        ContigSize_t cds_size = exon.second - exon.first;

        if (sub_sequence_offset < relative_offset + cds_size) {

          if (sub_sequence_offset <= relative_offset) {

            begin_offset = exon.first;

          } else {

            begin_offset = exon.first + (sub_sequence_offset - relative_offset);

          }

          if (sub_sequence_offset + sub_sequence_length > relative_offset + cds_size) {

            end_offset = exon.second;

          } else {

            end_offset = exon.second - ((relative_offset + cds_size) - (sub_sequence_offset + sub_sequence_length));

          }

          begin = base_sequence.getAlphabetString().begin() + (begin_offset - contig_offset);
          end = base_sequence.getAlphabetString().begin() + (end_offset - contig_offset);
          std::transform( begin, end, std::back_inserter(coding_sequence), convert_base);

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

      for (auto rit = exon_offset_map.rbegin(); rit != exon_offset_map.rend(); ++rit) {

        ContigSize_t cds_size = rit->second - rit->first;

        if (sub_sequence_offset < relative_offset + cds_size) {

          if (sub_sequence_offset <= relative_offset) {

            begin_offset = rit->second;

          } else {

            begin_offset = rit->second - (sub_sequence_offset - relative_offset);

          }

          if (sub_sequence_offset + sub_sequence_length > relative_offset + cds_size) {

            end_offset = rit->first;

          } else {

            end_offset = rit->first + ((relative_offset + cds_size) - (sub_sequence_offset + sub_sequence_length));

          }

          rbegin = base_sequence.getAlphabetString().rbegin() + (base_sequence.length() - (begin_offset - contig_offset));
          rend = base_sequence.getAlphabetString().rbegin() + (base_sequence.length() - (end_offset - contig_offset));
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

    ExecEnv::log().error("SequenceOffset::codingSubSequence(); Coding SubSequence length: {} NOT EQUAL to specified subsequence: {}",
                         coding_sequence.length(),
                         calculated_seq_size);

  }

  return DNA5SequenceCoding(std::move(coding_sequence), strand);

}

// Converts a DNA5SequenceLinear sequence to a DNA5SequenceCoding sequence
kgl::DNA5SequenceCoding kgl::SequenceOffset::codingSequence(const DNA5SequenceLinear& base_sequence, StrandSense strand) {


  StringCodingDNA5 coding_sequence;
  coding_sequence.reserve(base_sequence.length() + 1); // Just to make sure.


  switch(strand) {

    case StrandSense::FORWARD: {
      auto convert_base = [](DNA5::Alphabet base) { return DNA5::convertToCodingDNA5(base); };
      std::transform(base_sequence.getAlphabetString().begin(),
                     base_sequence.getAlphabetString().end(),
                     std::back_inserter(coding_sequence), convert_base);
    }
      break;

    case StrandSense::REVERSE: {
      auto complement_base = [](DNA5::Alphabet base) { return DNA5::complementNucleotide(base); };
      std::transform(base_sequence.getAlphabetString().rbegin(),
                     base_sequence.getAlphabetString().rend(),
                     std::back_inserter(coding_sequence), complement_base);
    }
      break;

  }

  return DNA5SequenceCoding(std::move(coding_sequence), strand);

}




// Returns bool false if contig_offset is not within the coding sequence defined by the coding_seq_ptr.
// If the contig_offset is in the coding sequence then a valid sequence_offset and the sequence length is returned.
// The offset is adjusted for strand type; the offset arithmetic is reversed for -ve strand sequences.
bool kgl::SequenceOffset::refOffsetWithinCodingSequence(const std::shared_ptr<const TranscriptionSequence>& coding_seq_ptr,
                                                        ContigOffset_t contig_offset,
                                                        ContigOffset_t &coding_sequence_offset,
                                                        ContigSize_t &coding_sequence_length) {

  StrandSense strand;
  ExonOffsetMap exon_offset_map;
  if (not exonOffsetAdapter(coding_seq_ptr, strand, exon_offset_map)) {

    ExecEnv::log().error("refOffsetWithinCodingSequence(), Contig: {}, Gene: {},  cannot generate EXON map for contig offset: {}",
                         coding_seq_ptr->contig()->contigId(), coding_seq_ptr->getGene()->id(), contig_offset);
    return false;
  }

  return offsetWithinCodingSequence(exon_offset_map, coding_seq_ptr->strand(), contig_offset, 0, coding_sequence_offset, coding_sequence_length);

}




// Returns bool false if contig_offset is not within the coding sequence defined by the coding_seq_ptr.
// If the contig_offset is in the coding sequence then a valid sequence_offset and the sequence length is returned.
// The offset is adjusted for strand type; the offset arithmetic is reversed for -ve strand sequences.
// Assumes exon maps are half intervals [begin,end)
bool kgl::SequenceOffset::offsetWithinCodingSequence(const ExonOffsetMap& exon_offset_map,
                                                     StrandSense strand,
                                                     ContigOffset_t sequence_offset,
                                                     ContigOffset_t sequence_contig_offset,
                                                     ContigOffset_t &coding_sequence_offset,
                                                     ContigSize_t &coding_sequence_length) {

  bool iscoding = false;
  ContigOffset_t coding_offset = 0;
  ContigSize_t coding_size = 0;

  if (exon_offset_map.empty()) {

    coding_sequence_offset = 0;
    coding_sequence_length = 0;
    return false;

  }

  ContigOffset_t contig_offset = sequence_offset + sequence_contig_offset;


  switch(strand) {

    case StrandSense::FORWARD: {

      for (auto exon : exon_offset_map) {

        // within the CDS
        if (contig_offset >= exon.first and contig_offset < exon.second) {

          coding_offset += (contig_offset - exon.first);
          iscoding = true;

        } else if (contig_offset >= exon.second) {

          coding_offset += (exon.second - exon.first);

        }

        coding_size += (exon.second - exon.first);

      }

    }
      break;

      // Careful with this logic as the CDS offsets are [begin, end). Therefore we must adjust the offset by -1.
    case StrandSense::REVERSE: {

      for (auto rit = exon_offset_map.rbegin(); rit != exon_offset_map.rend(); ++rit) {

        // within the CDS
        if (contig_offset >= rit->first and contig_offset < rit->second) {

          coding_offset += (rit->second - contig_offset) - 1;  // Careful here, remember [begin, end)
          iscoding = true;

        } else if (contig_offset < rit->first) {

          coding_offset += (rit->second - rit->first);

        }

        coding_size += (rit->second - rit->first);

      }

    }
      break;

  } // switch

  if (iscoding) {

    coding_sequence_offset = coding_offset;
    coding_sequence_length = coding_size;

  } else {

    coding_sequence_offset = 0;
    coding_sequence_length = 0;

  }

  return iscoding;

}


// Inverse of the above function. Given the stranded base offset within a coding sequence, return the corresponding contig offset.
// Returns bool false if then coding sequence offset is not within the coding sequence defined by the coding_seq_ptr.
// If the coding sequence_offset is in the coding sequence then a valid contig_offset and the sequence length is returned.
// The contig offset is adjusted for strand type; the offset arithmetic is reversed for -ve strand sequences.

bool kgl::SequenceOffset::refCodingSequenceContigOffset(const std::shared_ptr<const TranscriptionSequence>& coding_seq_ptr,
                                                        ContigOffset_t coding_sequence_offset,
                                                        ContigOffset_t &contig_offset,
                                                        ContigSize_t &coding_sequence_length) {

  StrandSense strand;
  ExonOffsetMap exon_offset_map;
  if (not exonOffsetAdapter(coding_seq_ptr, strand, exon_offset_map)) {

    ExecEnv::log().error("refCodingSequenceContigOffset(), Contig: {}, Gene: {},  cannot generate EXON map for CODING sequence offset: {}",
                         coding_seq_ptr->contig()->contigId(), coding_seq_ptr->getGene()->id(), contig_offset);
    return false;

  }

  return codingSequenceContigOffset(exon_offset_map, strand, coding_sequence_offset, contig_offset, coding_sequence_length);

}

// Returns bool false if the sequence_offset is not within the coding sequence defined by the coding_seq_ptr.
// If the sequence_offset is in the coding sequence then a valid contig_offset and the sequence length is returned.
// The contig offset is adjusted for strand type; the offset arithmetic is reversed for -ve strand sequences.
// Assumes exon maps are half intervals [begin,end)
bool kgl::SequenceOffset::codingSequenceContigOffset(const ExonOffsetMap& exon_offset_map,
                                                     StrandSense strand,
                                                     ContigOffset_t sequence_offset,
                                                     ContigOffset_t &contig_offset,
                                                     ContigSize_t &coding_sequence_length) {

  bool iscoding = false;
  ContigSize_t coding_size = 0;

  if (exon_offset_map.empty()) {

    contig_offset = 0;
    coding_sequence_length = 0;
    ExecEnv::log().error("SequenceOffset::codingSequenceContigOffset(), coding sequence with no exon regions");
    return false;

  }


  switch(strand) {

    case StrandSense::FORWARD: {

      ContigOffset_t begin_offset = 0;

      for (auto exon : exon_offset_map) {

        ContigOffset_t end_offset = (exon.second - exon.first) + begin_offset;
        // within the CDS
        if (sequence_offset >= begin_offset and sequence_offset < end_offset) {

          contig_offset = exon.first + sequence_offset - begin_offset;
          iscoding = true;

        } else if (sequence_offset >= end_offset) {

          begin_offset += (exon.second - exon.first);

        }

        coding_size += (exon.second - exon.first);

      }

    }
      break;

      // Careful with this logic as the CDS offsets are [begin, end). Therefore we must adjust the offset by -1.
    case StrandSense::REVERSE: {

      ContigOffset_t begin_offset = 0;

      for (auto rit = exon_offset_map.rbegin(); rit != exon_offset_map.rend(); ++rit) {

        ContigOffset_t end_offset = (rit->second - rit->first) + begin_offset;

        // within the CDS
        if (sequence_offset >= begin_offset and sequence_offset < end_offset) {

          contig_offset = (rit->second - 1) - (sequence_offset - begin_offset);  // Careful here, remember [begin, end)
          iscoding = true;

        } else if (sequence_offset >= end_offset) {

          begin_offset += (rit->second - rit->first);

        }

        coding_size += (rit->second - rit->first);

      }

    }
      break;

  } // switch

  if (iscoding) {

    coding_sequence_length = coding_size;

  } else {

    ExecEnv::log().error("SequenceOffset::codingSequenceContigOffset(), Sequence Offset: {} not in coding sequence", sequence_offset);
    contig_offset = 0;
    coding_sequence_length = 0;

  }

  return iscoding;

}

