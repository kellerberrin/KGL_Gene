//
// Created by kellerberrin on 1/12/18.
//

#ifndef KGL_GENOME_AUX_H
#define KGL_GENOME_AUX_H


#include <memory>
#include <string>
#include <vector>
#include <map>
#include "kgl_genome_types.h"
#include "kgl_sequence_amino.h"
#include "kgl_genome_feature.h"


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// These objects contain auxillary genome information sources e.g. TSS, Histone modification, Promoter sites etc.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


namespace kellerberrin {   //  organization level namespace
namespace genome {   // project level namespace




/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// A container to hold genome features.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class StructuredFeatures {

public:

  StructuredFeatures() = default;
  StructuredFeatures(const StructuredFeatures&) = default;
  virtual ~StructuredFeatures() = default;

  // Make this a virtual object.
  virtual const std::string featureType() const = 0;
  // Checks the feature type before adding.
  virtual bool checkAddFeature(std::shared_ptr<Feature>& feature_ptr) = 0;

  // False if not found.
  bool findFeatureId(const FeatureIdent_t& feature_id, std::vector<std::shared_ptr<Feature>>& feature_ptr_vec) const;

  const OffsetFeatureMap& offsetFeatureMap() const { return offset_feature_map_; }
  const IdFeatureMap& idFeatureMap() const { return id_feature_map_; }


protected:

  // Indexes by id and offset.
  void addFeature(std::shared_ptr<Feature>& feature_ptr);
  void verifyFeatureHierarchy();
  void clearHierarchy();

private:

  OffsetFeatureMap offset_feature_map_;
  IdFeatureMap id_feature_map_;

  void verifyContigOverlap();
  void verifySubFeatureDuplicates();
  void verifySuperFeatureDuplicates();
  void removeSubFeatureDuplicates();
  void removeSuperFeatureDuplicates();

};


using StructuredFeatureMap = std::map<std::string, std::shared_ptr<StructuredFeatures>>;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// A specialized container GENE and exon genome features.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class GeneExonFeatures : public StructuredFeatures {

public:

  GeneExonFeatures() = default;
  GeneExonFeatures(const GeneExonFeatures&) = default;
  ~GeneExonFeatures() override = default;

  const std::string featureType() const override { return GENE_EXON_FEATURE_; }

  void setupVerifyHierarchy();

  bool findGenes(ContigOffset_t offset, GeneVector &gene_ptr_vec) const;

  // Given a gene id and an mRNA (sequence id) return the CDS coding sequence.
  bool getCodingSequence(const FeatureIdent_t& gene_id,
                         const FeatureIdent_t& sequence_id,
                         std::shared_ptr<const CodingSequence>& coding_sequence_ptr) const;

  // Checks the feature type before adding (not TSS).
  bool checkAddFeature(std::shared_ptr<Feature>& feature_ptr) override;

  static constexpr const char* GENE_EXON_FEATURE_{"GeneExonFeatures"};


private:

  GeneMap gene_map_;

  void verifySubFeatureSuperFeatureDimensions();
  void createGeneMap();
  void setupFeatureHierarchy();
  void verifyGeneExonHierarchy();

};



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// A specialized container for Transcription Start Sequence (TSS) features.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class AdjalleyTSSFeatures : public StructuredFeatures {

public:

  AdjalleyTSSFeatures() = default;
  AdjalleyTSSFeatures(const AdjalleyTSSFeatures&) = default;
  ~AdjalleyTSSFeatures() override = default;

  const std::string featureType() const override { return ADJALLEY_TSS_FEATURE_; }

  void setupVerifyHierarchy(const StructuredFeatures& gene_super_features);

  // Checks the feature type before adding (must be TSS).
  bool checkAddFeature(std::shared_ptr<Feature>& feature_ptr) override;

  // Return all TSS features in this contig.
  TSSVector getTSSVector() const;

  static constexpr const char* ADJALLEY_TSS_FEATURE_{"AdjalleyTSSFeatures"};

private:

  void setupFeatureHierarchy(const StructuredFeatures& gene_super_features);


};






}   // namespace genome
}   // namespace kellerberrin


#endif //KGL_GENOME_AUX_H
