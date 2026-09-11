//
// Created by kellerberrin on 10/10/17.
//

#include "kel_exec_env.h"
#include "kgl_genome_feature.h"
#include "kgl_genome_contig.h"

#include <format>
#include <algorithm>
#include <ranges>

namespace kgl = kellerberrin::genome;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Feature members.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void kgl::Feature::addSubFeature(const FeatureIdent_t& sub_feature_id, std::shared_ptr<const Feature> sub_feature_ptr) {

  sub_features_.emplace(sub_feature_id, std::move(sub_feature_ptr));

}


void kgl::Feature::recursivePrintSubFeatures(size_t feature_level) const {

  ExecEnv::log().info("Level: {};  {}", feature_level, featureText());

  for (const auto& [feature_id, feature_ptr] : sub_features_) {

    feature_ptr->recursivePrintSubFeatures(feature_level + 1); // Recursive call for the sub-feature.

  }

}


std::string kgl::Feature::featureText(char delimiter) const {

  const std::string super_feature_id = hasSuperfeature() ? getSuperFeature()->id() : "<TopLevel>";

  return std::format(
      "Contig Id:{0}{1}{0}Feature Id:{0}{2}{0}Type:{0}{3}{0}SuperFeature:{0}{4}{0}SubFeatures:{0}{5}{0}Length:{0}{6}{0}Offset:[{0}{7}{0}{8}{0}) Strand:{0}{9}{0}Phase:{0}{10}{0}Description:{0}{11}",
      delimiter,
      contig_ref_ptr()->contigId(),
      id(),
      type(),
      super_feature_id,
      subFeatures().size(),
      sequence().length(),
      sequence().begin(),
      sequence().end(),
      sequence().strandText(),
      sequence().phase(),
      descriptionText(delimiter));

}

std::string kgl::Feature::descriptionText(char delimiter) const {

  std::string description_text;

  for (auto const& description : getAttributes().getDescription()) {

    description_text += description;
    description_text += delimiter;

  }

  return description_text;

}


bool kgl::Feature::equivalent(const Feature& lhs) const {

  return id_ == lhs.id_
         and type_ == lhs.type_
         and sequence_.equivalent(lhs.sequence_)
         and attributes_.equivalent(lhs.attributes_);

}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Gene Feature members.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::shared_ptr<const kgl::TranscriptionSequenceArray>
kgl::GeneFeature::getTranscriptionSequences(const std::shared_ptr<const GeneFeature>& gene_ptr) {

  auto sequence_array_ptr = std::make_shared<TranscriptionSequenceArray>();
  if (not getCodingSequences(gene_ptr, gene_ptr, *sequence_array_ptr)) {

    ExecEnv::log().error("GeneFeature::getTranscriptionSequences; Unable to retrieve coding sequences for Gene: {}, type: {}", gene_ptr->id(), gene_ptr->type());

  }

  return sequence_array_ptr;

}



// This routine is recursive. Assumes all the CDS and EXONs are leaf nodes in the feature hierarchy.
bool kgl::GeneFeature::getCodingSequences(const std::shared_ptr<const GeneFeature>& gene_ptr,
                                          const std::shared_ptr<const Feature>& parent_ptr,
                                          TranscriptionSequenceArray& sequence_array) {

  bool result = true;

  // Check if bottom level exists.
  if (gene_ptr->subFeatures().empty()) {

    return true;  // Allow empty genes.

  }

  std::vector<std::shared_ptr<const Feature>> leaf_features;
  bool CDS_found{false};
  // Loop through all sub_features.
  for (const auto& [feature_id, sub_feature_ptr] : parent_ptr->subFeatures()) {

    // If subfeature is bottom level.
    if (sub_feature_ptr->subFeatures().empty()) {

      if (sub_feature_ptr->superType() == CDS_TYPE_) {

        CDS_found = true;

      }
      leaf_features.push_back(sub_feature_ptr);

    } else { // Assume feature is a higher feature such as mRNA and recursively call this function for sub-features.

      result = result and getCodingSequences(gene_ptr, sub_feature_ptr, sequence_array);

    }

  } // for loop

  if (CDS_found) {
  // Remove all non-cds (EXONS) from the list
    std::erase_if(leaf_features, [](const std::shared_ptr<const Feature>& feature) { return feature->superType() != CDS_TYPE_; });

  } else {
  // Check that all the leaf features are EXONS and the gene is ncRNA.
    auto leaf_type = [](const std::shared_ptr<const Feature>& feature)->bool {
      return feature->superType() == EXON_TYPE_ or feature->superType() == ENHANCER_TYPE_;
    };

    if (not std::ranges::all_of(leaf_features, leaf_type)) {

      ExecEnv::log().warn("GeneFeature::getCodingSequences; Gene: {}, Unexpected leaf feature detected", gene_ptr->id());
      gene_ptr->recursivePrintSubFeatures();

    }

  }

  TranscriptionFeatureMap parent_feature_map;
  for (auto const& leaf : leaf_features) {

    auto [iter, insert_result] = parent_feature_map.emplace(leaf->sequence().begin(), leaf);

    // Some GFF files may have multiple coding features at the same logical level and the same begin offset.
    // This is true of the GFF supplied by NCBI for the SARS-COV-2 organism with multiple gene coding for the RdRp gene.
    if (not insert_result) {

      ExecEnv::log().warn("GeneFeature::getCodingSequences; Gene: {}, Duplicate coding feature: {}",
                          gene_ptr->id(),
                          leaf->featureText());
      gene_ptr->recursivePrintSubFeatures();
      return true;

    }

  }

  if (not parent_feature_map.empty()) {

    auto coding_sequence = std::make_shared<const TranscriptionSequence>(gene_ptr, parent_ptr, std::move(parent_feature_map));
    result = result and sequence_array.insertSequence(std::move(coding_sequence));

  }

  return result;

}

// Recursively find a feature super type in a feature hierarchy.
// If CDS is found then a protein gene, else ncRNA gene.
bool kgl::GeneFeature::findSuperType(const FeatureType_t& super_type, const std::shared_ptr<const Feature>& feature_ptr) {

  if (feature_ptr->superType() == super_type) {

    return true;

  }

  for (auto const& [feature_id, sub_feature_ptr] : feature_ptr->subFeatures()) {

    if (findSuperType(super_type, sub_feature_ptr)) {

      return true;

    }

  }

  return false;

}
