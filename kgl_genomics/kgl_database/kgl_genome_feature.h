///
// Created by kellerberrin on 10/10/17.
//

#ifndef KGL_GENOME_FEATURE_H
#define KGL_GENOME_FEATURE_H

#include "kgl_genome_attributes.h"
#include "kgl_genome_prelim.h"


namespace kellerberrin::genome {   //  organization level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Feature - Annotated contig features
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class ContigReference; // Forward decl;
using SubFeatureMap = std::multimap<const FeatureIdent_t, std::shared_ptr<const Feature>>;
using SuperFeaturePtr = std::weak_ptr<const Feature>;

class Feature {

public:

  Feature(FeatureIdent_t id,
          FeatureType_t type, // The actual Gff classification such as 'NCRNA_GENE'
          FeatureType_t super_type, // High level feature classification such as 'GENE'
          std::shared_ptr<const ContigReference> contig_ptr,
          FeatureSequence sequence): id_(std::move(id)), type_(std::move(type)), super_type_(std::move(super_type)),
                                     contig_ptr_(std::move(contig_ptr)), sequence_(sequence) {}
  Feature(const Feature&) = default;
  virtual ~Feature() = default;

//  Feature& operator=(const Feature&) = default;

  [[nodiscard]] const FeatureIdent_t& id() const { return id_; }
  [[nodiscard]] const FeatureSequence& sequence() const { return sequence_; }
  void sequence(const FeatureSequence& new_sequence) { sequence_ = new_sequence; }
  void setAttributes(const Attributes& attributes) { attributes_ = attributes; }
  [[nodiscard]] const Attributes& getAttributes() const { return attributes_; }
  [[nodiscard]] const FeatureType_t& type() const { return type_; } // The actual Gff classification such as 'NCRNA_GENE'
  [[nodiscard]] const FeatureType_t& superType() const { return super_type_; } // High level feature classification such as 'GENE'
  [[nodiscard]] bool verifyStrand(const TranscriptionFeatureMap& Feature_map) const;   // Check feature strand consistency
  [[nodiscard]] bool verifyCDSPhase(const TranscriptionSequenceArray& coding_seq_array) const; // Check the CDS phase for -ve and +ve strand genes
  void recusivelyPrintsubfeatures(size_t feature_level = 1) const; // useful debug function.
  [[nodiscard]] std::string featureText(char delimiter = ' ') const; // also useful for debug
  [[nodiscard]] std::string descriptionText(char delimiter = ' ') const; // Displays feature description.
  // Hierarchy routines.
  void clearHierachy() { sub_features_.clear(); super_feature_ptr_.reset(); }
  void addSubFeature(const FeatureIdent_t& sub_feature_id, std::shared_ptr<const Feature> sub_feature_ptr);

  constexpr static const char GENE_TYPE_[] = "GENE";
  constexpr static const char MRNA_TYPE_[] = "MRNA";
  constexpr static const char EXON_TYPE_[] = "EXON";
  constexpr static const char CDS_TYPE_[] = "CDS";
  constexpr static const char UTR5_TYPE_[] = "FIVE_PRIME_UTR";
  constexpr static const char UTR3_TYPE_[] = "THREE_PRIME_UTR";
  constexpr static const char TSS_TYPE_[] = "TSS_BLOCK";
  constexpr static const char ENHANCER_TYPE_[] = "ENHANCER";

  [[nodiscard]] bool isGene() const { return superType() == GENE_TYPE_; }
  [[nodiscard]] bool ismRNA() const { return superType() == MRNA_TYPE_; }
  [[nodiscard]] bool isUTR5() const { return superType() == UTR5_TYPE_; }
  [[nodiscard]] bool isUTR3() const { return superType() == UTR3_TYPE_; }
  [[nodiscard]] bool isTSS() const { return superType() == TSS_TYPE_; }

  // Always check for a NULL pointer with super feature.
  [[nodiscard]] bool hasSuperfeature() const { return getSuperFeature() != nullptr; }
  [[nodiscard]] std::shared_ptr<const Feature> getSuperFeature() const { return super_feature_ptr_.lock(); }
  void setSuperFeature(const std::shared_ptr<const Feature>& shared_super_feature) { super_feature_ptr_ = shared_super_feature; }

  [[nodiscard]] const SubFeatureMap& subFeatures() const { return sub_features_; }
  [[nodiscard]] SubFeatureMap& subFeatures() { return sub_features_; }

  [[nodiscard]] std::shared_ptr<const ContigReference> contig() const { return contig_ptr_; }

  // Compares features - primarily used for testing
  [[nodiscard]] bool equivalent(const Feature& lhs) const;

private:

  const FeatureIdent_t id_;
  const FeatureType_t type_;  // The actual Gff classification such as 'NCRNA_GENE'
  const FeatureType_t super_type_; // High level feature classification such as 'GENE'
  std::shared_ptr<const ContigReference> contig_ptr_;
  FeatureSequence sequence_;
  SubFeatureMap sub_features_;
  SuperFeaturePtr super_feature_ptr_;
  Attributes attributes_;

  [[nodiscard]] bool verifyMod3(const TranscriptionFeatureMap& feature_map) const;
};


using OffsetFeatureMap = std::multimap<ContigOffset_t, std::shared_ptr<Feature>>; // Contig features indexed by offset.
using IdFeatureMap = std::multimap<FeatureIdent_t, std::shared_ptr<Feature>>; // Contig features indexed by ident.


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// A Gene - can have multiple coding CDSFeatures/mRNAFeature sequences
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class GeneFeature : public Feature {

public:

  GeneFeature(const FeatureIdent_t &id,
              std::string type,
              const std::shared_ptr<const ContigReference> &contig_ptr,
              const FeatureSequence &sequence) : Feature(id, std::move(type), GENE_TYPE_, contig_ptr, sequence) {}

  GeneFeature(const GeneFeature &) = default;
  ~GeneFeature() override = default;

  // Public function returns a sequence map.
  [[nodiscard]] static std::unique_ptr<const TranscriptionSequenceArray> getTranscriptionSequences(const std::shared_ptr<const GeneFeature>& gene);

  constexpr static const char CODING_GENE_[] = "GENE";
  constexpr static const char PROTEIN_CODING_GENE_[] = "PROTEIN_CODING_GENE";
  constexpr static const char NCRNA_GENE_[] = "NCRNA_GENE";

  [[nodiscard]] static bool proteinCoding(const std::shared_ptr<const GeneFeature>& gene_ptr) { return findSuperType(CDS_TYPE_, gene_ptr); }
  [[nodiscard]] static bool ncRNACoding(const std::shared_ptr<const GeneFeature>& gene_ptr) { return not proteinCoding(gene_ptr); }

private:

  // Initializes the map_ptr with a list of coding sequences found for the gene. Recursive.
  // Specifying an optional identifier is used
  [[nodiscard]] static bool getCodingSequences(const std::shared_ptr<const GeneFeature>& gene,
                                               const std::shared_ptr<const Feature>& parent_ptr,
                                               TranscriptionSequenceArray& sequence_array_ptr);

  // Recursively find a feature super type in a feature hierarchy.
  // If CDS is found then a protein gene, else ncRNA gene.
  [[nodiscard]] static bool findSuperType(const FeatureType_t& type, const std::shared_ptr<const Feature>& feature);

};

using GeneVector = std::vector<std::shared_ptr<const GeneFeature>>;  // Multiple alternative genes for sequence region.
using GeneMap = std::multimap<ContigOffset_t, std::shared_ptr<const GeneFeature>>;  // Inserted using the END offset as key.





}   // end namespace


#endif //KGL_GENOME_FEATURE_H
