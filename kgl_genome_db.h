//
// Created by kellerberrin on 7/10/17.
//

#ifndef KGL_GENOME_DB_H
#define KGL_GENOME_DB_H


#include <memory>
#include <string>
#include <vector>
#include <map>
#include "kgl_genome_types.h"
#include "kgl_sequence_amino.h"
#include "kgl_genome_feature.h"
#include "kgl_mt_data.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ContigFeatures - A contiguous region, the associated sequence, and all features that map onto that region/sequence.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


using OffsetFeatureMap = std::multimap<ContigOffset_t, std::shared_ptr<Feature>>; // Contig features indexed by offset.
using IdFeatureMap = std::multimap<FeatureIdent_t, std::shared_ptr<Feature>>; // Contig features indexed by ident.
using CDSArray = std::vector<std::shared_ptr<CDSFeature>>;
using GeneMap = std::multimap<ContigOffset_t, std::shared_ptr<GeneFeature>>;  // Inserted using the END offset as key.
using GeneVector = std::vector<std::shared_ptr<GeneFeature>>;  // Multiple alternative genes for sequence region.


class ContigFeatures {

public:

  ContigFeatures(const ContigId_t &contig_id,
                 std::shared_ptr<DNA5Sequence> sequence_ptr) : contig_id_(contig_id), sequence_ptr_(sequence_ptr) {}
  ContigFeatures(const ContigFeatures&) = default;
  ~ContigFeatures() = default;

  ContigFeatures& operator=(const ContigFeatures&) = default;

  bool addFeature(std::shared_ptr<Feature>& feature_ptr);
  // false if not found.
  bool findFeatureId(FeatureIdent_t& feature_id, std::vector<std::shared_ptr<Feature>>& feature_ptr_vec);
  // false if offset is not in an cds else returns a vector of cds (these will be in different genes).
  bool findOffsetCDS(ContigOffset_t offset, CDSArray& cds_array) const;
  // false if offset is not in a gene, else (true) returns a vector of ptrs to the genes.
  bool findGenes(ContigOffset_t offset, GeneVector &gene_ptr_vec) const;

  bool setTranslationTable(size_t table) { return coding_sequence_.settranslationTable(table); }
  std::string translationTableName() const { return coding_sequence_.translationTableName(); }

  const ContigId_t& contigId() const { return contig_id_; }
  const DNA5Sequence& sequence() const { return *sequence_ptr_; }
  ContigSize_t contigSize() { return sequence_ptr_->length(); }

  void setupFeatureHierarchy();
  void verifyFeatureHierarchy();
  void verifyCDSPhasePeptide();

private:

  ContigId_t contig_id_;
  std::shared_ptr<DNA5Sequence> sequence_ptr_;
  OffsetFeatureMap offset_feature_map_;
  IdFeatureMap id_feature_map_;
  GeneMap gene_map_;
  CodingSequenceDNA5 coding_sequence_;  // Amino Acid translation table, can be unique for contig (e.g. mitochondria)

  void verifyContigOverlap();
  void verifySubFeatureSuperFeatureDimensions();
  void verifySubFeatureDuplicates();
  void verifySuperFeatureDuplicates();
  void removeSubFeatureDuplicates();
  void removeSuperFeatureDuplicates();
  void createGeneMap();
  // Check all gene coding sequences for start and end codons and nonsense (intermediate stop codon) mutations.
  bool verifyCodingSequences(const SortedCDSVector& sorted_cds_vec) const;

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// GenomeDatabase - A map of contigs
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


using GenomeSequenceMap = std::map<ContigId_t, std::shared_ptr<ContigFeatures>>;
class GenomeDatabase {

public:

  explicit GenomeDatabase() {}
  GenomeDatabase(const GenomeDatabase&) = default;
  ~GenomeDatabase() = default;

  GenomeDatabase& operator=(const GenomeDatabase&) = default;

  // Return false if contig already exists.
  bool addContigSequence(const ContigId_t& contig, std::shared_ptr<DNA5Sequence> sequence_ptr);
  // Returns false if key not found.
  bool getContigSequence(const ContigId_t& contig, std::shared_ptr<ContigFeatures>& contig_ptr) const;

  void createVerifyGenomeDatabase();

  void setTranslationTable(size_t table);

  void registerContigData(std::shared_ptr<ContigCountData>& contig_data_ptr) const;

private:

  GenomeSequenceMap genome_sequence_map_;

  void setupFeatureHierarchy();
  void verifyFeatureHierarchy();

};


}   // namespace genome
}   // namespace kellerberrin


#endif //KGL_GENOME_DB_H
