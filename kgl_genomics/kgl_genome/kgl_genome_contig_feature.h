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



namespace kellerberrin::genome {   //  organization::project level namespace



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
  [[nodiscard]] virtual std::string featureType() const = 0;
  // Checks the feature type before adding.
  [[nodiscard]] virtual bool checkAddFeature(std::shared_ptr<Feature>& feature_ptr) = 0;

  // False if not found.
  [[nodiscard]] std::vector<std::shared_ptr<const Feature>> findFeatureId(const FeatureIdent_t& feature_id) const;

  [[nodiscard]] const OffsetFeatureMap& offsetFeatureMap() const { return offset_feature_map_; }
  [[nodiscard]] const IdFeatureMap& idFeatureMap() const { return id_feature_map_; }

  // Compares contig feature structures - primarily used for testing
  [[nodiscard]] bool equivalent(const StructuredFeatures& lhs) const;

protected:

  // Indexes by id and offset.
  void addFeature(std::shared_ptr<Feature>& feature_ptr);
  void verifyFeatureHierarchy();
  void clearHierarchy();

private:

  OffsetFeatureMap offset_feature_map_;
  IdFeatureMap id_feature_map_;

  void verifyContigOverlap() const;
  size_t verifySubFeatureDuplicates() const;
  void removeSubFeatureDuplicates();

};


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// A specialized container for GENE and exon genome features.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class GeneExonFeatures : public StructuredFeatures {

public:

  GeneExonFeatures() = default;
  GeneExonFeatures(const GeneExonFeatures&) = default;
  ~GeneExonFeatures() override = default;

  [[nodiscard]] std::string featureType() const override { return GENE_EXON_FEATURE_; }

  void setupVerifyHierarchy();

  [[nodiscard]] const GeneMap& geneMap() const { return gene_map_; }

  // Checks the feature type before adding.
  [[nodiscard]] bool checkAddFeature(std::shared_ptr<Feature>& feature_ptr) override;

  static constexpr const char* GENE_EXON_FEATURE_{"GeneExonFeatures"};


private:

  GeneMap gene_map_;

  void verifySubFeatureSuperFeatureDimensions();
  void createGeneMap();
  void setupFeatureHierarchy();
  void verifyGeneExonHierarchy();
  static bool checkSubFeatures( const std::shared_ptr<const Feature>& feature_ptr
                              , const std::shared_ptr<const Feature>& sub_feature_ptr);
  static bool checkSuperFeature( const std::shared_ptr<const Feature>& feature_ptr
                               , const std::shared_ptr<const Feature>& super_feature_ptr);

  };




}   // end namespace


#endif //KGL_GENOME_AUX_H
