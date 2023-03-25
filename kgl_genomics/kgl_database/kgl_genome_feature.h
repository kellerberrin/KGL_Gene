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
          FeatureType_t type,
          std::shared_ptr<const ContigReference> contig_ptr,
          FeatureSequence sequence): id_(std::move(id)), type_(std::move(type)), contig_ptr_(std::move(contig_ptr)), sequence_(sequence) {}
  Feature(const Feature&) = default;
  virtual ~Feature() = default;

//  Feature& operator=(const Feature&) = default;

  [[nodiscard]] const FeatureIdent_t& id() const { return id_; }
  [[nodiscard]] const FeatureSequence& sequence() const { return sequence_; }
  void sequence(const FeatureSequence& new_sequence) { sequence_ = new_sequence; }
  void setAttributes(const Attributes& attributes) { attributes_ = attributes; }
  [[nodiscard]] const Attributes& getAttributes() const { return attributes_; }
  [[nodiscard]] const FeatureType_t& featureType() const { return type_; }
  [[nodiscard]] bool verifyStrand(const SortedCDS& sorted_cds) const;   // Check feature strand consistency
  [[nodiscard]] bool verifyCDSPhase(std::shared_ptr<const CodingSequenceArray> coding_seq_ptr) const; // Check the CDS phase for -ve and +ve strand genes
  [[nodiscard]] std::shared_ptr<const Feature> getGene() const; // returns null pointer if not found
  void recusivelyPrintsubfeatures(long feature_level = 1) const; // useful debug function.
  [[nodiscard]] std::string featureText() const; // also useful for debug
  // Hierarchy routines.
  void clearHierachy() { sub_features_.clear(); super_feature_ptr_.reset(); }
  void addSubFeature(const FeatureIdent_t& sub_feature_id, std::shared_ptr<const Feature> sub_feature_ptr);

  constexpr static const char GENE_TYPE_[] = "GENE";
  constexpr static const char MRNA_TYPE_[] = "MRNA";
  constexpr static const char CDS_TYPE[] = "CDS";
  constexpr static const char UTR5_TYPE[] = "FIVE_PRIME_UTR";
  constexpr static const char UTR3_TYPE[] = "THREE_PRIME_UTR";
  constexpr static const char TSS_TYPE[] = "TSS_BLOCK";

  [[nodiscard]] bool isGene() const { return type_ == GENE_TYPE_; }
  [[nodiscard]] bool ismRNA() const { return type_ == MRNA_TYPE_; }
  [[nodiscard]] bool isCDS() const { return type_ == CDS_TYPE; }
  [[nodiscard]] bool isUTR5() const { return type_ == UTR5_TYPE; }
  [[nodiscard]] bool isUTR3() const { return type_ == UTR3_TYPE; }
  [[nodiscard]] bool isTSS() const { return type_ == TSS_TYPE; }

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
  const FeatureType_t type_;
  std::shared_ptr<const ContigReference> contig_ptr_;
  FeatureSequence sequence_;
  SubFeatureMap sub_features_;
  SuperFeaturePtr super_feature_ptr_;
  Attributes attributes_;

  [[nodiscard]] static bool verifyMod3(const SortedCDS& sorted_cds);
};


using OffsetFeatureMap = std::multimap<ContigOffset_t, std::shared_ptr<Feature>>; // Contig features indexed by offset.
using IdFeatureMap = std::multimap<FeatureIdent_t, std::shared_ptr<Feature>>; // Contig features indexed by ident.


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// A Gene - can have multiple coding CDSFeatures/mRNAFeature sequences
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class GeneFeature : public Feature {

public:

  GeneFeature(const FeatureIdent_t &id,
              const std::shared_ptr<const ContigReference> &contig_ptr,
              const FeatureSequence &sequence) : Feature(id, GENE_TYPE_, contig_ptr, sequence) {}

  GeneFeature(const GeneFeature &) = default;
  ~GeneFeature() override = default;

  // Public function returns a sequence map.
  [[nodiscard]] static std::shared_ptr<const CodingSequenceArray> getCodingSequences(std::shared_ptr<const GeneFeature> gene);

private:

  // Initializes the map_ptr with a list of coding sequences found for the gene. Recursive.
  // Specifying an optional identifier is used
  [[nodiscard]] static bool getCodingSequences( std::shared_ptr<const GeneFeature> gene,
                                                std::shared_ptr<const Feature> cds_parent,
                                                std::shared_ptr<CodingSequenceArray>& sequence_array_ptr,
                                                const FeatureIdent_t& optional_identifier = FeatureIdent_t(""));


};

using GeneVector = std::vector<std::shared_ptr<const GeneFeature>>;  // Multiple alternative genes for sequence region.
using GeneMap = std::multimap<ContigOffset_t, std::shared_ptr<const GeneFeature>>;  // Inserted using the END offset as key.



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// A Coding CDS Feature
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class CDSFeature : public Feature {

public:

  CDSFeature(const FeatureIdent_t &id,
             const CDSPhaseType_t phase,
             const std::shared_ptr<const ContigReference> &contig_ptr,
             const FeatureSequence &sequence) : Feature(id, CDS_TYPE, contig_ptr, sequence), phase_(phase) {}

  CDSFeature(const CDSFeature &) = default;
  ~CDSFeature() override = default;

  [[nodiscard]] CDSPhaseType_t CDSphase() const { return phase_; }

  constexpr static const char* CDS_TYPE = "CDS";

private:

  CDSPhaseType_t phase_;

};



}   // end namespace


#endif //KGL_GENOME_FEATURE_H
