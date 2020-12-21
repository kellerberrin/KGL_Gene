//
// Created by kellerberrin on 3/12/18.
//

#ifndef KGL_GENOME_CONTIG_AUX_H
#define KGL_GENOME_CONTIG_AUX_H

#include "kgl_genome_contig_feature.h"


namespace kellerberrin::genome {   //  organization::project level namespace


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

  [[nodiscard]] const std::string featureType() const override { return ADJALLEY_TSS_FEATURE_; }

  void setupVerifyHierarchy(const StructuredFeatures& gene_super_features);

  // Checks the feature type before adding (must be TSS).
  [[nodiscard]] bool checkAddFeature(std::shared_ptr<Feature>& feature_ptr) override;

  // Return all TSS features in this contig.
  [[nodiscard]] TSSVector getTSSVector() const;

  static constexpr const char* ADJALLEY_TSS_FEATURE_{"AdjalleyTSSFeatures"};

private:

  void setupFeatureHierarchy(const StructuredFeatures& gene_super_features);


};



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// A collection of objects that contain auxillary genome information sources e.g. TSS, Histone modification, Promoter sites etc.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


class AuxContigFeatures {

public:

  AuxContigFeatures() = default;
  AuxContigFeatures(const AuxContigFeatures&) = default;
  ~AuxContigFeatures() = default;

  const AdjalleyTSSFeatures& getTSSfeatures() const { return adjalley_TSS_Features_; }

  void setupVerifyHierarchy(const StructuredFeatures& gene_super_features);

  void checkAddFeature(std::shared_ptr<Feature>& feature_ptr);

private:

  AdjalleyTSSFeatures adjalley_TSS_Features_;

};



}   // end namespace



#endif //KGL_GENOME_CONTIG_AUX_H
