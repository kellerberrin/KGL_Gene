///
// Created by kellerberrin on 10/10/17.
//

#ifndef KGL_GENOME_FEATURE_H
#define KGL_GENOME_FEATURE_H

#include "kgl_attributes.h"
#include "kgl_genome_prelim.h"


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Feature - Annotated contig features
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class ContigFeatures; // Forward decl;
using SubFeatureMap = std::multimap<const FeatureIdent_t, std::shared_ptr<Feature>>;
using SuperFeatureMap = std::multimap<const FeatureIdent_t, std::shared_ptr<Feature>>;

class Feature {

public:

  Feature(const FeatureIdent_t& id,
          const FeatureType_t& type,
          const std::shared_ptr<const ContigFeatures>& contig_ptr,
          const FeatureSequence& sequence): id_(id), type_(type), contig_ptr_(contig_ptr), sequence_(sequence) {}
  Feature(const Feature&) = default;
  virtual ~Feature() = default;

  Feature& operator=(const Feature&) = default;

  const FeatureIdent_t& id() const { return id_; }
  const FeatureSequence& sequence() const { return sequence_; }
  void sequence(const FeatureSequence& new_sequence) { sequence_ = new_sequence; }
  void setAttributes(const Attributes& attributes) { attributes_ = attributes; }
  const Attributes& getAttributes() const { return attributes_; }
  const FeatureType_t& featureType() const { return type_; }
  bool verifyStrand(const SortedCDS& sorted_cds) const;   // Check feature strand consistency
  bool verifyCDSPhase(std::shared_ptr<const CodingSequenceArray> coding_seq_ptr) const; // Check the CDS phase for -ve and +ve strand genes
  std::shared_ptr<Feature> getGene() const; // returns null pointer if not found
  void recusivelyPrintsubfeatures(long feature_level = 1) const; // useful debug function.
  // Hierarchy routines.
  void clearHierachy() { sub_features_.clear(); super_features_.clear(); }
  void addSuperFeature(const FeatureIdent_t &super_feature_id, const std::shared_ptr<Feature> &super_feature_ptr);
  void addSubFeature(const FeatureIdent_t& sub_feature_id, const std::shared_ptr<Feature>& sub_feature_ptr);
  virtual bool isCDS() const { return false; }
  virtual bool isGene() const { return false; }
  virtual bool ismRNA() const { return false; }
  virtual bool isEXON() const { return false; }

  const SuperFeatureMap& superFeatures() const { return super_features_; }
  const SubFeatureMap& subFeatures() const { return sub_features_; }
  SuperFeatureMap& superFeatures() { return super_features_; }
  SubFeatureMap& subFeatures() { return sub_features_; }

  std::shared_ptr<const ContigFeatures> contig() const { return contig_ptr_; }

private:

  FeatureIdent_t id_;
  FeatureType_t type_;
  std::shared_ptr<const ContigFeatures> contig_ptr_;
  FeatureSequence sequence_;
  SubFeatureMap sub_features_;
  SuperFeatureMap super_features_;
  Attributes attributes_;

  bool verifyMod3(const SortedCDS& sorted_cds) const;
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// A Coding CDS Feature
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class CDSFeature : public Feature {

public:

  CDSFeature(const FeatureIdent_t &id,
            const CDSPhaseType_t phase,
            const std::shared_ptr<const ContigFeatures> &contig_ptr,
            const FeatureSequence &sequence) : Feature(id, CDS_TYPE, contig_ptr, sequence), phase_(phase) {}

  CDSFeature(const CDSFeature &) = default;
  ~CDSFeature() override = default;

  CDSFeature &operator=(const CDSFeature &) = default;

  CDSPhaseType_t phase() const { return phase_; }
  void phase(CDSPhaseType_t _phase) { phase_ = _phase; }

  virtual bool isCDS() const final { return true; }
  // CDS Type.
  constexpr static const char* CDS_TYPE = "CDS";

private:

  CDSPhaseType_t phase_;

};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// A EXONFeature - co-occurs with the CDS features.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class EXONFeature : public Feature {

public:

  EXONFeature(const FeatureIdent_t &id,
              const std::shared_ptr<const ContigFeatures> &contig_ptr,
              const FeatureSequence &sequence) : Feature(id, EXON_TYPE, contig_ptr, sequence) {}

  EXONFeature(const EXONFeature &) = default;
  ~EXONFeature() override = default;

  EXONFeature &operator=(const EXONFeature &) = default;

  bool isEXON() const final { return true; }

  // EXON Type.
  constexpr static const char* EXON_TYPE = "EXON";

private:


};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// A mRNAFeature - generally the parent of the coding CDS Features
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class mRNAFeature : public Feature {

public:

  mRNAFeature(const FeatureIdent_t &id,
              const std::shared_ptr<const ContigFeatures> &contig_ptr,
              const FeatureSequence &sequence) : Feature(id, MRNA_TYPE, contig_ptr, sequence) {}

  mRNAFeature(const mRNAFeature &) = default;
  ~mRNAFeature() override = default;

  mRNAFeature &operator=(const mRNAFeature &) = default;

  bool ismRNA() const final { return true; }

  // MRNA Type.
  constexpr static const char* MRNA_TYPE = "MRNA";

private:


};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// A Gene - can have multiple coding CDSFeatures/mRNAFeature sequences
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

using GeneVector = std::vector<std::shared_ptr<const GeneFeature>>;  // Multiple alternative genes for sequence region.
class GeneFeature : public Feature {

public:

  GeneFeature(const FeatureIdent_t &id,
             const std::shared_ptr<const ContigFeatures> &contig_ptr,
             const FeatureSequence &sequence) : Feature(id, GENE_TYPE, contig_ptr, sequence) {}

  GeneFeature(const GeneFeature &) = default;
  ~GeneFeature() override = default;

  GeneFeature &operator=(const GeneFeature &) = default;

  // Public function returns a sequence map.
  static std::shared_ptr<const CodingSequenceArray> getCodingSequences(std::shared_ptr<const GeneFeature> gene);

  bool isGene() const final { return true; }

  // GENE Type.
  constexpr static const char* GENE_TYPE = "GENE";

private:

  // Initializes the map_ptr with a list of coding sequences found for the gene. Recursive.
  static bool getCodingSequences(std::shared_ptr<const GeneFeature> gene,
                                 std::shared_ptr<const Feature> cds_parent,
                                 std::shared_ptr<CodingSequenceArray>& sequence_array_ptr);


};


}   // namespace genome
}   // namespace kellerberrin


#endif //KGL_GENOME_FEATURE_H
