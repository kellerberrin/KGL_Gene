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
using GeneMap = std::multimap<ContigOffset_t, std::shared_ptr<GeneFeature>>;  // Inserted using the END offset as key.


class ContigFeatures {

public:

  ContigFeatures(const ContigId_t& contig_id,
                 std::shared_ptr<DNA5SequenceContig> sequence_ptr) : contig_id_(contig_id),
                                                                     sequence_ptr_(sequence_ptr) {}
  ContigFeatures(const ContigFeatures&) = default;
  ~ContigFeatures() = default;

  ContigFeatures& operator=(const ContigFeatures&) = default;

  bool addFeature(std::shared_ptr<Feature>& feature_ptr);
  // false if not found.
  bool findFeatureId(const FeatureIdent_t& feature_id, std::vector<std::shared_ptr<Feature>>& feature_ptr_vec) const;
  // false if offset is not in a gene, else (true) returns a vector of ptrs to the genes.
  bool findGenes(ContigOffset_t offset, GeneVector &gene_ptr_vec) const;
  const GeneMap& getGeneMap() const { return gene_map_; }

  bool setTranslationTable(size_t table) { return coding_sequence_.settranslationTable(table); }
  std::string translationTableName() const { return coding_sequence_.translationTableName(); }

  const ContigId_t& contigId() const { return contig_id_; }
  const DNA5SequenceContig& sequence() const { return *sequence_ptr_; }
  ContigSize_t contigSize() { return sequence_ptr_->length(); }

  void setupFeatureHierarchy();
  void verifyFeatureHierarchy();
  void verifyCDSPhasePeptide();

  // Given a gene id and an mRNA (sequence id) return the CDS coding sequence.
  bool getCodingSequence(const FeatureIdent_t& gene_id,
                         const FeatureIdent_t& sequence_id,
                         std::shared_ptr<const CodingSequence>& coding_sequence_ptr) const;

  // Given a CDS coding sequence, return the corresponding DNA base sequence (strand adjusted).
  bool getDNA5SequenceCoding(const std::shared_ptr<const CodingSequence>& coding_sequence_ptr,
                             std::shared_ptr<DNA5SequenceCoding>& sequence_ptr) const;

  // Generate Amino acid sequences using the table specified for this contig.
  std::shared_ptr<AminoSequence> getAminoSequence(std::shared_ptr<const CodingSequence> coding_seq_ptr) const;
  std::shared_ptr<AminoSequence> getAminoSequence(std::shared_ptr<const DNA5SequenceCoding> sequence_ptr) const;
  AminoAcid::Alphabet getAminoAcid(const Codon& codon) const { return coding_sequence_.getAmino(codon); }


private:

  ContigId_t contig_id_;
  std::shared_ptr<DNA5SequenceContig> sequence_ptr_;
  OffsetFeatureMap offset_feature_map_;
  IdFeatureMap id_feature_map_;
  GeneMap gene_map_;
  TranslateToAmino coding_sequence_;  // Amino Acid translation table, unique for contig (e.g. mitochondria)

  void verifyContigOverlap();
  void verifySubFeatureSuperFeatureDimensions();
  void verifySubFeatureDuplicates();
  void verifySuperFeatureDuplicates();
  void removeSubFeatureDuplicates();
  void removeSuperFeatureDuplicates();
  void createGeneMap();
  // Check all gene coding sequences for start and end codons and nonsense (intermediate stop codon) mutations.
  bool verifyCodingSequences(const std::shared_ptr<const GeneFeature> gene_ptr,
                             std::shared_ptr<const CodingSequenceArray> coding_seq_ptr) const;

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
  bool addContigSequence(const ContigId_t& contig, std::shared_ptr<DNA5SequenceContig> sequence_ptr);
  // Returns false if key not found.
  bool getContigSequence(const ContigId_t& contig, std::shared_ptr<ContigFeatures>& contig_ptr) const;

  void createVerifyGenomeDatabase();

  void setTranslationTable(size_t table);

  void registerContigData(std::shared_ptr<ContigCountData>& contig_data_ptr) const;

  const GenomeSequenceMap& getMap() const { return genome_sequence_map_; }

  size_t contigCount() const { return getMap().size(); }

private:

  GenomeSequenceMap genome_sequence_map_;

  void setupFeatureHierarchy();
  void verifyFeatureHierarchy();

};


}   // namespace genome
}   // namespace kellerberrin


#endif //KGL_GENOME_DB_H
