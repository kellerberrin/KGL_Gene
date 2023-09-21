//
// Created by kellerberrin on 29/05/23.
//

#ifndef KGL_GENOME_INTERVAL_H
#define KGL_GENOME_INTERVAL_H


#include "kgl_genome_genome.h"
#include "kgl_variant_db.h"
#include "kel_interval_unsigned.h"
#include "kel_interval_set.h"
#include "kel_interval_map.h"



namespace kellerberrin::genome {   //  organization::project level namespace




/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// A gene interval node.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// A map of named gene coding transcripts.
using GeneCodingTranscriptMap = std::map<std::string, IntervalSetLower>;
class GeneIntervalStructure {

public:

  explicit GeneIntervalStructure(const std::shared_ptr<const GeneFeature>& gene_feature)
  : gene_interval_(gene_feature->sequence().begin(), gene_feature->sequence().end()) { codingInterval(gene_feature); }
  ~GeneIntervalStructure() = default;

  GeneIntervalStructure(const GeneIntervalStructure& copy) = default;

  // Object access
  [[nodiscard]] const std::shared_ptr<const GeneFeature>& getGene() const { return gene_feature_; }
  [[nodiscard]] const GeneCodingTranscriptMap& codingTranscripts() const { return gene_coding_transcripts_; }
  [[nodiscard]] const OpenRightUnsigned& geneInterval() const { return gene_interval_; }
  [[nodiscard]] const IntervalSetLower& transcriptUnion() const { return transcript_union_; }

  // Given an offset, does the offset fall within a defined gene interval (can include 5 prime, coding intervals, introns, 3 prime).
  [[nodiscard]] bool isWithinGene(const Variant& variant) const { return gene_interval_.containsOffset(variant.offset()); }
  // Test if a variant modifies any of the CODING transcripts of this gene structure.
  [[nodiscard]] bool codingModifier(const Variant& variant) const;
  // Test if a variant modifies the specified transcript.
  [[nodiscard]] bool transcriptModifier(const Variant& variant, const FeatureIdent_t& Transcript) const;

private:

  std::shared_ptr<const GeneFeature> gene_feature_;
  GeneCodingTranscriptMap gene_coding_transcripts_;
  IntervalSetLower transcript_union_;
  OpenRightUnsigned gene_interval_;

  void codingInterval(const std::shared_ptr<const GeneFeature>& gene_vector);
  // Used with the functions above to determine if a contig + offset resides within a gene interval or the coding intervals of a gene.
  [[nodiscard]] bool isSameContig(const ContigId_t& contig) const { return contig == (gene_feature_->contig()->contigId()); }

};


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Implemented using the interval multimap as gene intervals can (and do) overlap.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


using ContigIntervalMap = std::map<ContigId_t, IntervalLowerMultiMap <std::shared_ptr<const GeneIntervalStructure>>>;
class IntervalCodingVariants {

public:

  explicit IntervalCodingVariants(const std::vector<std::shared_ptr<const GeneFeature>>& gene_vector) { InitializeGeneVector(gene_vector); }
  explicit IntervalCodingVariants(const std::shared_ptr<const GenomeReference>& reference_ptr);
  ~IntervalCodingVariants() = default;

  // Returns true if the variant is within a gene coding region.
  [[nodiscard]] bool codingRegionVariant(const Variant& variant) const;

  [[nodiscard]] const ContigIntervalMap& getCodingMap() const { return contig_interval_map_; }

private:

  ContigIntervalMap contig_interval_map_;

  void InitializeGeneVector(const std::vector<std::shared_ptr<const GeneFeature>>& gene_vector);

};



} // Namespace



#endif //KGL_GENOME_INTERVAL_H
