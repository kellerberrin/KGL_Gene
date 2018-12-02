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

  explicit StructuredFeatures(const std::string& feature_type) : feature_type_(feature_type) {}
  StructuredFeatures(const StructuredFeatures&) = default;
  virtual ~StructuredFeatures() = default;

  const std::string& featureType() const { return feature_type_; }

  bool addFeature(std::shared_ptr<Feature>& feature_ptr);

  // false if not found.
  bool findFeatureId(const FeatureIdent_t& feature_id, std::vector<std::shared_ptr<Feature>>& feature_ptr_vec) const;

  const OffsetFeatureMap offetFeatureMap() const { return offset_feature_map_; }

private:

  const std::string feature_type_;
  OffsetFeatureMap offset_feature_map_;
  IdFeatureMap id_feature_map_;

};


using StruturedFeatureMap = std::map<std::string, std::shared_ptr<StructuredFeatures>>;



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// A specialized container GENE and exon genome features.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class GeneExonFeatures : public StructuredFeatures {

public:

  GeneExonFeatures() : StructuredFeatures(GENE_EXON_FEATURE_) {}
  GeneExonFeatures(const GeneExonFeatures&) = default;
  ~GeneExonFeatures() = default;

  static constexpr const char* GENE_EXON_FEATURE_{"GeneExonFeatures"};

private:


};



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// A specialized container for Transcription Start Sequence (TSS) features.
//
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class TSSFeatures : public StructuredFeatures {

public:

  TSSFeatures() : StructuredFeatures(TSS_FEATURE_) {}
  TSSFeatures(const TSSFeatures&) = default;
  ~TSSFeatures() = default;

  static constexpr const char* TSS_FEATURE_{"TSSFeatures"};

  // Return all TSS features in this contig.
  TSSVector getTSSVector() const;

private:


};






}   // namespace genome
}   // namespace kellerberrin


#endif //KGL_GENOME_AUX_H
