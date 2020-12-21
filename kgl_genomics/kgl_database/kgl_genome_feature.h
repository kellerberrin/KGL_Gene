///
// Created by kellerberrin on 10/10/17.
//

#ifndef KGL_GENOME_FEATURE_H
#define KGL_GENOME_FEATURE_H

#include "kgl_database/kgl_attributes.h"
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

  Feature(const FeatureIdent_t& id,
          const FeatureType_t& type,
          const std::shared_ptr<const ContigReference>& contig_ptr,
          const FeatureSequence& sequence): id_(id), type_(type), contig_ptr_(contig_ptr), sequence_(sequence) {}
  Feature(const Feature&) = default;
  virtual ~Feature() = default;

  Feature& operator=(const Feature&) = default;

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
  [[nodiscard]] virtual bool isCDS() const { return false; }
  [[nodiscard]] virtual bool isGene() const { return false; }
  [[nodiscard]] virtual bool ismRNA() const { return false; }
  [[nodiscard]] virtual bool isPSEUDOGENE() const { return false; }
  [[nodiscard]] virtual bool isTSS() const { return false; }
  [[nodiscard]] virtual bool isEXON() const { return false; }


  // Always check for a NULL pointer with super feature.
  [[nodiscard]] bool hasSuperfeature() const { return getSuperFeature() != nullptr; }
  [[nodiscard]] std::shared_ptr<const Feature> getSuperFeature() const { return super_feature_ptr_.lock(); }
  void setSuperFeature(std::shared_ptr<const Feature> shared_super_feature) { super_feature_ptr_ = shared_super_feature; }

  [[nodiscard]] const SubFeatureMap& subFeatures() const { return sub_features_; }
  [[nodiscard]] SubFeatureMap& subFeatures() { return sub_features_; }

  [[nodiscard]] std::shared_ptr<const ContigReference> contig() const { return contig_ptr_; }

private:

  FeatureIdent_t id_;
  FeatureType_t type_;
  std::shared_ptr<const ContigReference> contig_ptr_;
  FeatureSequence sequence_;
  SubFeatureMap sub_features_;
  SuperFeaturePtr super_feature_ptr_;
  Attributes attributes_;

  [[nodiscard]] bool verifyMod3(const SortedCDS& sorted_cds) const;
};


using OffsetFeatureMap = std::multimap<ContigOffset_t, std::shared_ptr<Feature>>; // Contig features indexed by offset.
using IdFeatureMap = std::multimap<FeatureIdent_t, std::shared_ptr<Feature>>; // Contig features indexed by ident.


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

  CDSFeature &operator=(const CDSFeature &) = default;

  [[nodiscard]] CDSPhaseType_t phase() const { return phase_; }
  void phase(CDSPhaseType_t _phase) { phase_ = _phase; }

  [[nodiscard]] virtual bool isCDS() const final { return true; }
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
              const std::shared_ptr<const ContigReference> &contig_ptr,
              const FeatureSequence &sequence) : Feature(id, EXON_TYPE, contig_ptr, sequence) {}

  EXONFeature(const EXONFeature &) = default;
  ~EXONFeature() override = default;

  EXONFeature &operator=(const EXONFeature &) = default;

  [[nodiscard]] bool isEXON() const final { return true; }

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
              const std::shared_ptr<const ContigReference> &contig_ptr,
              const FeatureSequence &sequence) : Feature(id, MRNA_TYPE, contig_ptr, sequence) {}

  mRNAFeature(const mRNAFeature &) = default;
  ~mRNAFeature() override = default;

  mRNAFeature &operator=(const mRNAFeature &) = default;

  [[nodiscard]] bool ismRNA() const final { return true; }

  // MRNA Type.
  constexpr static const char* MRNA_TYPE = "MRNA";

private:


};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// A Pseudogene feature may or may not have CDS or EXON subfeatures.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class PSEUDOGENEFeature : public Feature {

public:

  PSEUDOGENEFeature(const FeatureIdent_t &id,
                    const std::shared_ptr<const ContigReference> &contig_ptr,
                    const FeatureSequence &sequence) : Feature(id, PSEUDOGENE_TYPE, contig_ptr, sequence) {}

  PSEUDOGENEFeature(const PSEUDOGENEFeature &) = default;
  ~PSEUDOGENEFeature() override = default;

  PSEUDOGENEFeature &operator=(const PSEUDOGENEFeature &) = default;

  [[nodiscard]] bool isPSEUDOGENE() const final { return true; }

  // EXON Type.
  constexpr static const char* PSEUDOGENE_TYPE = "PSEUDOGENE";

private:


};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// A Gene - can have multiple coding CDSFeatures/mRNAFeature sequences
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class GeneFeature : public Feature {

public:

  GeneFeature(const FeatureIdent_t &id,
             const std::shared_ptr<const ContigReference> &contig_ptr,
             const FeatureSequence &sequence) : Feature(id, GENE_TYPE, contig_ptr, sequence) {}

  GeneFeature(const GeneFeature &) = default;
  ~GeneFeature() override = default;

  GeneFeature &operator=(const GeneFeature &) = default;

  // Public function returns a sequence map.
  [[nodiscard]] static std::shared_ptr<const CodingSequenceArray> getCodingSequences(std::shared_ptr<const GeneFeature> gene);

  [[nodiscard]] bool isGene() const final { return true; }

  // GENE Type.
  constexpr static const char* GENE_TYPE = "GENE";

private:

  // Initializes the map_ptr with a list of coding sequences found for the gene. Recursive.
  // Specifying an optional identifier is used 
  [[nodiscard]] static bool getCodingSequences( std::shared_ptr<const GeneFeature> gene,
                                                std::shared_ptr<const Feature> cds_parent,
                                                std::shared_ptr<CodingSequenceArray>& sequence_array_ptr,
                                                const FeatureIdent_t& optional_identifier = FeatureIdent_t(""));


};

using GeneVector = std::vector<std::shared_ptr<const GeneFeature>>;  // Multiple alternative genes for sequence region.
using GeneMap = std::multimap<ContigOffset_t, std::shared_ptr<GeneFeature>>;  // Inserted using the END offset as key.


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// A Transcription Start Sequence, can be multiple TSS per gene.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

class TSSFeature;
using TSSVector = std::vector<std::shared_ptr<const TSSFeature>>;

class TSSFeature : public Feature {

public:

  TSSFeature(const FeatureIdent_t &id,
             const std::shared_ptr<const ContigReference> &contig_ptr,
             const FeatureSequence &sequence) : Feature(id, TSS_TYPE, contig_ptr, sequence) {}

  TSSFeature(const TSSFeature &) = default;
  ~TSSFeature() override = default;

  TSSFeature &operator=(const TSSFeature &) = default;

  [[nodiscard]] bool isTSS() const final { return true; }

  // TSS Type.
  constexpr static const char* TSS_TYPE = "TSS_BLOCK";
  // Unassigned to a feature.
  constexpr static const char* TSS_UNASSIGNED = "NewTranscript";

private:


};




}   // end namespace


#endif //KGL_GENOME_FEATURE_H
