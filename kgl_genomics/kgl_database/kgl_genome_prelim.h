//
// Created by kellerberrin on 7/11/17.
//

#ifndef KGL_GENOME_PRELIM_H
#define KGL_GENOME_PRELIM_H

#include <memory>
#include <string>
#include <vector>
#include <map>
#include "kgl_genome_types.h"


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
  // They always correspond to zero based offsets on the relevant contig.
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

enum class TranscriptionSequenceType { PROTEIN, NCRNA, EMPTY};
// Check a Protein transcription sequence for validity. Note that all ncRNA sequences are trivially valid.
enum class ProteinSequenceValidity { VALID, EMPTY,  NOT_MOD3, NO_START_CODON, NONSENSE_MUTATION, NO_STOP_CODON };
using TranscriptionFeatureMap = std::map<ContigOffset_t, std::shared_ptr<const Feature>>;
class TranscriptionSequence {

public:

  TranscriptionSequence(std::shared_ptr<const GeneFeature> gene_ptr,
                        std::shared_ptr<const Feature> parent_ptr,
                        TranscriptionFeatureMap feature_map): gene_ptr_(std::move(gene_ptr)),
                                                              parent_ptr_(std::move(parent_ptr)),
                                                              transcription_feature_map_(std::move(feature_map)) {}
  ~TranscriptionSequence() = default;

  [[nodiscard]] const TranscriptionFeatureMap& getFeatureMap() const { return transcription_feature_map_; }
  [[nodiscard]] size_t codingFeatures() const { return getFeatureMap().size(); }
  [[nodiscard]] std::shared_ptr<const ContigReference> contig() const;
  [[nodiscard]] std::shared_ptr<const GeneFeature> getGene() const { return gene_ptr_; }
  [[nodiscard]] std::shared_ptr<const Feature> getParent() const { return parent_ptr_; }
  [[nodiscard]] StrandSense strand() const;
  [[nodiscard]] ContigOffset_t prime_5() const; // Offset of the 5 prime nucleotide closest to the sequence (strand adjusted start - 1).
  void prime_5_region(ContigSize_t requested_size, ContigOffset_t& begin_offset, ContigSize_t& size) const;
  [[nodiscard]] ContigOffset_t prime_3() const; // Offset of the 3 prime nucleotide closest to the sequence (strand adjusted end).
  void prime_3_region(ContigSize_t requested_size, ContigOffset_t& begin_offset, ContigSize_t& size) const;
  [[nodiscard]] ContigOffset_t start() const; // Offset of the start of the sequence - not strand adjusted.
  [[nodiscard]] ContigOffset_t end() const; // Offset of the end of the sequence (last nucleotide + 1) - not strand adjusted.
  [[nodiscard]] ContigSize_t codingNucleotides() const; // Total number of nucleotides in all CDS.
  [[nodiscard]] TranscriptionSequenceType codingType() const;
  [[nodiscard]] static ProteinSequenceValidity checkValidProtein( const std::shared_ptr<const TranscriptionSequence>& coding_seq_ptr,
                                                                  bool verbose = false); // ncRNA are trivially valid.

private:

  std::shared_ptr<const GeneFeature> gene_ptr_;
  std::shared_ptr<const Feature> parent_ptr_;  // Generally an mRNAFeature, or whatever was the CDS superFeature().
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
  std::unique_ptr<const TranscriptionSequenceArray> validTranscriptionArray() const;

private:

  TranscriptionSequenceMap transcription_sequence_map_;

};



}   // end namespace


#endif //KGL_GENOME_PRELIM_H
