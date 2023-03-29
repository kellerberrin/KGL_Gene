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
// CodingSequence - The Gene and mRNA feature (if present) and the CDS features necessary for a protein sequence.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class Feature; // Forward decl.
class GeneFeature; // Forward decl.
class ContigReference; // Forward decl.

enum class TranscriptionSequenceType { PROTEIN, NCRNA, EMPTY};
using TranscribedFeatureMap = std::map<ContigOffset_t, std::shared_ptr<const Feature>>;
class CodingSequence {

public:

  CodingSequence(std::shared_ptr<const GeneFeature> gene_ptr,
                 std::shared_ptr<const Feature> parent_ptr,
                 TranscribedFeatureMap feature_map): gene_ptr_(std::move(gene_ptr)),
                                                     parent_ptr_(std::move(parent_ptr)),
                                                     coding_feature_map_(std::move(feature_map)) {}
  ~CodingSequence() = default;

  [[nodiscard]] const TranscribedFeatureMap& getFeatureMap() const { return coding_feature_map_; }
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

private:

  std::shared_ptr<const GeneFeature> gene_ptr_;
  std::shared_ptr<const Feature> parent_ptr_;  // Generally an mRNAFeature, or whatever was the CDS superFeature().
  TranscribedFeatureMap coding_feature_map_;

};

// A sorted array of coding sequences. Sorted by CDS parent (generally mRNA) feature ident.
// #define CODING_SEQUENCE_ISMAP 1  // Uncomment this if the requirement is for distinct CDS parent features.
#ifdef CODING_SEQUENCE_ISMAP
using TranscriptionSequenceMap = std::map<const FeatureIdent_t, std::shared_ptr<const CodingSequence>>; // For distinct parent features.
#else
using TranscriptionSequenceMap = std::multimap<const FeatureIdent_t, std::shared_ptr<const CodingSequence>>;  // Same parent features permitted.
#endif

class CodingSequenceArray {

public:

  CodingSequenceArray() = default;
  ~CodingSequenceArray() = default;

  [[nodiscard]] const TranscriptionSequenceMap& getMap() const { return coding_sequence_map_; }
  [[nodiscard]] TranscriptionSequenceMap& getMap() { return coding_sequence_map_; }
  [[nodiscard]] TranscriptionSequenceType codingType() const; // Assumes that all coding sequences are the same type.

  [[nodiscard]] bool insertCodingSequence(std::shared_ptr<const CodingSequence> coding_sequence);

  [[nodiscard]] size_t size() const { return coding_sequence_map_.size(); }
  [[nodiscard]] bool empty() const { return size() == 0; }

  [[nodiscard]] std::shared_ptr<const CodingSequence> getFirst() const;

  static void printSequence(std::shared_ptr<const CodingSequenceArray> coding_seq_ptr);

private:

  TranscriptionSequenceMap coding_sequence_map_;

};



}   // end namespace


#endif //KGL_GENOME_PRELIM_H
