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


enum class StrandSense { FORWARD = '+', REVERSE = '-'};
class FeatureSequence {

public:

  FeatureSequence(const ContigOffset_t& begin_offset,
                  const ContigOffset_t& end_offset,
                  const StrandSense& strand_sense)
  : begin_offset_(begin_offset),
    end_offset_(end_offset),
    strand_sense_(strand_sense) {}
  ~FeatureSequence() = default;
  FeatureSequence(const FeatureSequence&) = default;

  FeatureSequence& operator=(const FeatureSequence&) = default;

  ContigOffset_t begin() const { return begin_offset_; }
  ContigOffset_t end() const { return end_offset_; }
  StrandSense strand() const { return strand_sense_; }
  void begin(ContigOffset_t begin) { begin_offset_ = begin; }
  void end(ContigOffset_t end) { end_offset_ = end; }
  void strand(StrandSense strand) { strand_sense_ = strand; }


private:

  // Important - the offsets of the features use the half open interval [begin, end) convention and are ZERO based.
  // This is different from the gff format which are 1-based and are specified as the closed interval [begin, end].
  // Thus begin_offset_ = gff.begin - 1
  // And end_offset_ = gff.end
  // Note that always begin_offset_ < end_offset_. They are not strand adjusted.
  // They always correspond to zero based offsets on the relevant contig.
  ContigOffset_t begin_offset_;    // begin is the first nucleotide of the feature.
  ContigOffset_t end_offset_;      // end points past the last nucleotide of the feature; end = last_offset+1
  StrandSense strand_sense_;      // can be 'UNKNOWN'

};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// CodingSequence - The Gene and mRNA feature (if present) and the CDS features necessary for a protein sequence.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class Feature; // Forward decl.
class TSSFeature; // Forward decl
class CDSFeature; // Forward decl
class EXONFeature; // Forward decl
class GeneFeature; // Forward decl.
class mRNAFeature; // Forward decl.
class ContigFeatures; // Forward decl.

using SortedCDS = std::map<ContigOffset_t, std::shared_ptr<const CDSFeature>>;
class CodingSequence {

public:

  CodingSequence(std::shared_ptr<const GeneFeature> gene_ptr,
                 std::shared_ptr<const Feature> cds_parent_ptr,
                 SortedCDS sorted_cds): gene_ptr_(std::move(gene_ptr)),
                                        cds_parent_ptr_(std::move(cds_parent_ptr)),
                                        sorted_cds_(std::move(sorted_cds)) {}
  ~CodingSequence() = default;

  const SortedCDS& getSortedCDS() const { return sorted_cds_; }
  std::shared_ptr<const ContigFeatures> contig() const;
  std::shared_ptr<const GeneFeature> getGene() const { return gene_ptr_; }
  std::shared_ptr<const Feature> getCDSParent() const { return cds_parent_ptr_; }
  StrandSense strand() const;
  ContigOffset_t prime_5() const; // Offset of the 5 prime nucleotide closest to the sequence (strand adjusted start - 1).
  void prime_5_region(ContigSize_t requested_size, ContigOffset_t& begin_offset, ContigSize_t& size) const;
  ContigOffset_t prime_3() const; // Offset of the 3 prime nucleotide closest to the sequence (strand adjusted end).
  void prime_3_region(ContigSize_t requested_size, ContigOffset_t& begin_offset, ContigSize_t& size) const;
  ContigOffset_t start() const; // Offset of the start of the sequence - not strand adjusted.
  ContigOffset_t end() const; // Offset of the end of the sequence (last nucleotide + 1) - not strand adjusted.
  ContigSize_t codingNucleotides() const;

private:

  std::shared_ptr<const GeneFeature> gene_ptr_;
  std::shared_ptr<const Feature> cds_parent_ptr_;  // Generally an mRNAFeature, or whatever was the CDS superFeature().
  SortedCDS sorted_cds_;

};

// A sorted array of coding sequences. Sorted by CDS parent (mRNA) feature ident.
using CodingSequenceMap = std::map<const FeatureIdent_t, std::shared_ptr<const CodingSequence>>;
class CodingSequenceArray {

public:

  CodingSequenceArray() = default;
  ~CodingSequenceArray() = default;

  const CodingSequenceMap& getMap() const { return coding_sequence_map_; }
  CodingSequenceMap& getMap() { return coding_sequence_map_; }

  bool insertCodingSequence(std::shared_ptr<const CodingSequence> coding_sequence);

  size_t size() const { return coding_sequence_map_.size(); }
  bool empty() const { return size() == 0; }

  std::shared_ptr<const CodingSequence> getFirst() const;

  static void printCodingSequence(std::shared_ptr<const CodingSequenceArray> coding_seq_ptr);

private:

  CodingSequenceMap coding_sequence_map_;

};



}   // end namespace


#endif //KGL_GENOME_PRELIM_H
