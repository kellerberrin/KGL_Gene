//
// Created by kellerberrin on 7/11/17.
//

#ifndef KGL_GENOME_PRELIM_H
#define KGL_GENOME_PRELIM_H

#include "kgl_genome_types.h"

#include "kel_interval_unsigned.h"
#include "kel_interval_set.h"

#include <memory>
#include <string>
#include <vector>
#include <map>


namespace kellerberrin::genome {   //  organization level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// FeatureSequence - Feature location and strand (if applicable)
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

enum class StrandSense: char { FORWARD = '+', REVERSE = '-'};
class FeatureSequence {

public:

  FeatureSequence(const ContigOffset_t& begin_offset,
                  const ContigOffset_t& end_offset,
                  const StrandSense& strand_sense,
                  uint32_t phase = 0)
  : begin_offset_(begin_offset),
    end_offset_(end_offset),
    strand_sense_(strand_sense),
    phase_(phase) {}
  ~FeatureSequence() = default;
  FeatureSequence(const FeatureSequence&) = default;

  FeatureSequence& operator=(const FeatureSequence&) = default;

  [[nodiscard]] ContigOffset_t begin() const { return begin_offset_; }
  [[nodiscard]] ContigOffset_t end() const { return end_offset_; }
  [[nodiscard]] OpenRightUnsigned interval() const { return OpenRightUnsigned(begin(), end()); }
  [[nodiscard]] StrandSense strand() const { return strand_sense_; }
  [[nodiscard]] uint32_t phase() const { return phase_; }
  [[nodiscard]] ContigSize_t length() const { return end_offset_ - begin_offset_; }
  [[nodiscard]] char strandText() const { return static_cast<char>(strand()); }
  void begin(ContigOffset_t begin) { begin_offset_ = begin; }
  void end(ContigOffset_t end) { end_offset_ = end; }
  void strand(StrandSense strand) { strand_sense_ = strand; }
  // Returns the strand adjusted feature begin (-ve is end-1) to the
  // target strand adjusted feature begin (-ve is end-1)
  // and returns the relative begin transcription distance as a +ve offset.
  [[nodiscard]] ContigOffset_t distance(const FeatureSequence& compare_feature) const;

  // Primarily used for testing.
  [[nodiscard]] bool equivalent(const FeatureSequence& lhs) const { return begin() == lhs.begin() and end() == lhs.end() and strand() == lhs.strand(); }

private:

  // Important - the offsets of the features use the half open interval [begin, end) convention and are ZERO based.
  // This is different from the gff format which are 1-based and are specified as the closed interval [begin, end].
  // Thus begin_offset_ = gff.begin - 1
  // And end_offset_ = gff.end
  // Note that always begin_offset_ < end_offset_. They are not strand adjusted.
  // They always correspond to zero based offsets on the relevant contig_ref_ptr.
  ContigOffset_t begin_offset_;    // begin is the first nucleotide of the feature.
  ContigOffset_t end_offset_;      // end points past the last nucleotide of the feature; end = last_offset+1
  StrandSense strand_sense_;      // default is '+'
  uint32_t phase_;

};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// TranscriptionSequence - The Gene and mRNA feature (if present) and the CDS features necessary for a protein sequence.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class Feature; // Forward decl.
class GeneFeature; // Forward decl.
class ContigReference; // Forward decl.

// Check a transcription coding sequence for validity. Note that all ncRNA sequences are trivially valid.
enum class CodingSequenceValidity { NCRNA, VALID_PROTEIN, EMPTY,  NOT_MOD3, NO_START_CODON, NONSENSE_MUTATION, NO_STOP_CODON };

enum class TranscriptionSequenceType { PROTEIN, NCRNA, EMPTY};

using TranscriptionFeatureMap = std::map<ContigOffset_t, std::shared_ptr<const Feature>>;
class TranscriptionSequence {

public:

  TranscriptionSequence(std::shared_ptr<const GeneFeature> gene_ptr,
                        std::shared_ptr<const Feature> parent_ptr,
                        TranscriptionFeatureMap feature_map): gene_ptr_(std::move(gene_ptr)),
                                                              parent_ptr_(std::move(parent_ptr)),
                                                              transcription_feature_map_(std::move(feature_map)) {}
  ~TranscriptionSequence() = default;

  [[nodiscard]] IntervalSetLower getExonIntervals() const;
  [[nodiscard]] IntervalSetLower getIntronIntervals() const; // Empty for a 1 exon gene.
  [[nodiscard]] const TranscriptionFeatureMap& getFeatureMap() const { return transcription_feature_map_; }
  [[nodiscard]] std::shared_ptr<const ContigReference> contig() const;
  [[nodiscard]] std::shared_ptr<const GeneFeature> getGene() const { return gene_ptr_; }
  [[nodiscard]] std::shared_ptr<const Feature> getParent() const { return parent_ptr_; }
  [[nodiscard]] StrandSense strand() const;
  [[nodiscard]] OpenRightUnsigned prime5Region(ContigSize_t requested_size) const;
  [[nodiscard]] OpenRightUnsigned prime3Region(ContigSize_t requested_size) const;
  [[nodiscard]] ContigOffset_t start() const; // Zero-based offset [start, end) of the start of the sequence - not strand adjusted.
  [[nodiscard]] ContigOffset_t end() const; // Zero-based offset [start, end) of the end of the sequence (last nucleotide + 1) - not strand adjusted.
  [[nodiscard]] OpenRightUnsigned interval() const { return OpenRightUnsigned(start(), end()); }
  [[nodiscard]] ContigSize_t codingNucleotides() const; // Total number of nucleotides in all CDS.
  [[nodiscard]] TranscriptionSequenceType codingType() const;
  [[nodiscard]] static CodingSequenceValidity checkSequenceStatus(const std::shared_ptr<const TranscriptionSequence>& transcript_ptr);
  [[nodiscard]] static bool checkValidProtein(CodingSequenceValidity sequence_status) { return sequence_status == CodingSequenceValidity::VALID_PROTEIN; }
  // ncRNA sequences are trivially valid.
  [[nodiscard]] static bool checkValidSequence(CodingSequenceValidity sequence_status) {
    return checkValidProtein(sequence_status) or sequence_status == CodingSequenceValidity::NCRNA;
  }

private:

  std::shared_ptr<const GeneFeature> gene_ptr_;
  std::shared_ptr<const Feature> parent_ptr_;  // Generally an mRNAFeature, or whatever was the CDS array superFeature().
  TranscriptionFeatureMap transcription_feature_map_;

};

// A sorted array of coding sequences. Sorted by CDS parent (generally mRNA) feature ident.
// #define CODING_SEQUENCE_ISMAP 1  // Uncomment this if the requirement is for distinct CDS parent features.
#ifdef CODING_SEQUENCE_ISMAP
using TranscriptionSequenceMap = std::map<const FeatureIdent_t, std::shared_ptr<const TranscriptionSequence>>; // For distinct parent features.
#else
using TranscriptionSequenceMap = std::multimap<const FeatureIdent_t, std::shared_ptr<const TranscriptionSequence>>;  // Same parent features permitted.
#endif

class TranscriptionSequenceArray {

public:

  TranscriptionSequenceArray() = default;
  ~TranscriptionSequenceArray() = default;

  [[nodiscard]] const TranscriptionSequenceMap& getMap() const { return transcription_sequence_map_; }
  [[nodiscard]] TranscriptionSequenceMap& getMap() { return transcription_sequence_map_; }
  [[nodiscard]] TranscriptionSequenceType codingType() const; // Assumes that all coding sequences are the same type.

  [[nodiscard]] bool insertSequence(std::shared_ptr<const TranscriptionSequence> sequence_ptr);

  [[nodiscard]] size_t size() const { return transcription_sequence_map_.size(); }
  [[nodiscard]] bool empty() const { return transcription_sequence_map_.empty(); }

  [[nodiscard]] std::shared_ptr<const TranscriptionSequence> getFirst() const;

  static void printSequence(std::shared_ptr<const TranscriptionSequenceArray> sequence_ptr);

  // Only ncRNA and valid Protein transcriptions.
  [[nodiscard]] std::unique_ptr<const TranscriptionSequenceArray> validTranscriptionArray() const;

private:

  TranscriptionSequenceMap transcription_sequence_map_;

};

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
  [[nodiscard]] const std::vector<double>& nonsenseMutationSize() const { return nonsense_mutation_size_; }
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
  std::vector<double> nonsense_mutation_size_;
  size_t no_stop_codon_{0};


};

}   // end namespace


#endif //KGL_GENOME_PRELIM_H
