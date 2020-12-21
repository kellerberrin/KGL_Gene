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
#include "kgl_sequence/kgl_sequence_amino.h"
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
  [[nodiscard]] virtual const std::string featureType() const = 0;
  // Checks the feature type before adding.
  [[nodiscard]] virtual bool checkAddFeature(std::shared_ptr<Feature>& feature_ptr) = 0;

  // False if not found.
  [[nodiscard]] bool findFeatureId(const FeatureIdent_t& feature_id, std::vector<std::shared_ptr<const Feature>>& feature_ptr_vec) const;

  [[nodiscard]] const OffsetFeatureMap& offsetFeatureMap() const { return offset_feature_map_; }
  [[nodiscard]] const IdFeatureMap& idFeatureMap() const { return id_feature_map_; }


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
  void removeSubFeatureDuplicates();

};


using StructuredFeatureMap = std::map<std::string, std::shared_ptr<StructuredFeatures>>;


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

  [[nodiscard]] const std::string featureType() const override { return GENE_EXON_FEATURE_; }

  void setupVerifyHierarchy();

  [[nodiscard]] const GeneMap& geneMap() const { return gene_map_; }

  [[nodiscard]] bool findGenes(ContigOffset_t offset, GeneVector &gene_ptr_vec) const;

  // Given a gene id and an mRNA (sequence id) return the CDS coding sequence.
  [[nodiscard]] bool getCodingSequence( const FeatureIdent_t& gene_id,
                                        const FeatureIdent_t& sequence_id,
                                        std::shared_ptr<const CodingSequence>& coding_sequence_ptr) const;

  // Checks the feature type before adding (not TSS).
  [[nodiscard]] bool checkAddFeature(std::shared_ptr<Feature>& feature_ptr) override;

  static constexpr const char* GENE_EXON_FEATURE_{"GeneExonFeatures"};


private:

  GeneMap gene_map_;

  void verifySubFeatureSuperFeatureDimensions();
  void createGeneMap();
  void setupFeatureHierarchy();
  void verifyGeneExonHierarchy();

};




}   // end namespace


#endif //KGL_GENOME_AUX_H
