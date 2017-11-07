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


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// FeatureAttributes Object to hold { key=value } pairs.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


using AttributeMap = std::multimap<const std::string, std::string>;
class FeatureAttributes {

public:

  explicit FeatureAttributes() = default;
  FeatureAttributes(const FeatureAttributes&) = default;
  ~FeatureAttributes() = default;

  FeatureAttributes& operator=(const FeatureAttributes&) = default;

  // General access routines.
  bool getAttributes(const std::string& key, std::vector<std::string>& value_vec) const; // false if no key.
  void insertAttribute(const std::string& key, const std::string& value); // Always succeeds; keys are uppercase.
  void getAllAttributes(std::vector<std::pair<std::string, std::string>>& all_key_value_pairs) const;

  // Attribute keys.
  constexpr static const char* ID_KEY = "ID";
  constexpr static const char* SUPER_FEATURE_KEY = "PARENT";

  // Convenience access routines.
  bool getIds(std::vector<std::string> &value_vec) const { return getAttributes(ID_KEY, value_vec); }
  bool getSuperFeatureIds(std::vector<std::string> &value_vec) const { return getAttributes(SUPER_FEATURE_KEY, value_vec); }


private:

  AttributeMap attributes_;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// FeatureSequence - Feature location and strand (if applicable)
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


enum class StrandSense { FORWARD = '+', REVERSE = '-', UNKNOWN = '.'};
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

  ContigOffset_t begin_offset_;    // [O, size-1] based contig sequence offset.
  ContigOffset_t end_offset_;      // [O, size-1] based contig sequence offset.
  StrandSense strand_sense_;

};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// CodingSequence - The Gene and mRNA feature (if present) and the CDS features necessary for a protein sequence.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class Feature; // Forward decl.
class CDSFeature; // Forward decl
class EXONFeature; // Forward decl
class GeneFeature; // Forward decl.
class mRNAFeature; // Forward decl.

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
  std::shared_ptr<const GeneFeature> getGene() const { return gene_ptr_; }
  std::shared_ptr<const Feature> getCDSParent() const { return cds_parent_ptr_; }
  bool isWithinCoding(ContigOffset_t contig_offset) const;

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

  void mergeArrays(std::shared_ptr<const CodingSequenceArray> merge_array);


  static void printCodingSequence(std::shared_ptr<const CodingSequenceArray> coding_seq_ptr);

private:

  CodingSequenceMap coding_sequence_map_;

};



}   // namespace genome
}   // namespace kellerberrin




#endif //KGL_GENOME_PRELIM_H
