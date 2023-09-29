//
// Created by kellerberrin on 10/10/17.
//

#include "kel_exec_env.h"
#include "kgl_genome_feature.h"
#include "kgl_genome_contig.h"

#include <sstream>
#include <list>

namespace kgl = kellerberrin::genome;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Feature members.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void kgl::Feature::addSubFeature(const FeatureIdent_t& sub_feature_id, std::shared_ptr<const Feature> sub_feature_ptr) {

  sub_features_.insert(std::make_pair(sub_feature_id, sub_feature_ptr));

}


void kgl::Feature::recusivelyPrintsubfeatures(size_t feature_level) const {

  ExecEnv::log().info("Level: {};  {}", feature_level, featureText());

  for (const auto& [feature_id, feature_ptr] : sub_features_) {

    feature_ptr->recusivelyPrintsubfeatures(feature_level + 1); // Recursive call for the sub-feature.

  }

}


std::string kgl::Feature::featureText(char delimiter) const {

  std::stringstream ss;

  ss << "Contig Id:" << delimiter
     << contig_ref_ptr()->contigId() << delimiter
     << "Feature Id:" << delimiter
     << id() << delimiter
     << "Type:" << delimiter
     << type() << delimiter
     << "SuperFeature:" << delimiter
     << (hasSuperfeature() ? getSuperFeature()->id() : "<TopLevel>") << delimiter
     << "SubFeatures:" << delimiter
     << subFeatures().size() << delimiter
     << "Length:" << delimiter
     << sequence().length() << delimiter
     << "Offset:[" << delimiter
     << sequence().begin() << delimiter
     << sequence().end() << delimiter
     << ") Strand:" << delimiter
      << sequence().strandText() << delimiter
     << "Phase:" << delimiter
     << sequence().phase() << delimiter
     << "Description:" << delimiter
     << descriptionText(delimiter);

  return ss.str();

}

std::string kgl::Feature::descriptionText(char delimiter) const {

  std::string description_text;

  auto description_vector = getAttributes().getDescription();
  for (auto const& description : description_vector) {

    description_text += description;
    description_text += delimiter;

  } // Gene

  return description_text;

}


bool kgl::Feature::equivalent(const Feature& lhs) const {

  if (id_ != lhs.id_) {

    return false;

  }

  if (type_ != lhs.type_) {

    return false;

  }

  if (sequence_.equivalent(lhs.sequence_)) {

    return false;

  }

  if (attributes_.equivalent(lhs.attributes_)) {

    return false;

  }

  return true;

}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Gene Feature members.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::unique_ptr<const kgl::TranscriptionSequenceArray>
kgl::GeneFeature::getTranscriptionSequences(const std::shared_ptr<const GeneFeature>& gene_ptr) {

  std::unique_ptr<TranscriptionSequenceArray> sequence_array_ptr(std::make_unique<TranscriptionSequenceArray>());
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

  std::list<std::shared_ptr<const Feature>> leaf_features;
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
    leaf_features.remove_if([](const std::shared_ptr<const Feature>& feature)->bool { return feature->superType() != CDS_TYPE_; });

  } else {
  // Check that all the leaf features are EXONS and the gene is ncRNA.
    auto leaf_type = [](const std::shared_ptr<const Feature>& feature)->bool {
      return feature->superType() == EXON_TYPE_ or feature->superType() == ENHANCER_TYPE_;
    };

    if (not std::ranges::all_of(leaf_features, leaf_type)) {

      ExecEnv::log().warn("GeneFeature::getCodingSequences; Gene: {}, Unexpected leaf feature detected", gene_ptr->id());
      gene_ptr->recusivelyPrintsubfeatures();

    }

  }

  TranscriptionFeatureMap parent_feature_map;
  if (not leaf_features.empty()) {

    for (auto const& leaf : leaf_features) {

      auto [iter, insert_result] = parent_feature_map.insert(std::make_pair(leaf->sequence().begin(), leaf));

      // Some GFF files may have multiple coding features at the same logical level and the same begin offset.
      // This is true of the GFF supplied by NCBI for the SARS-COV-2 organism with multiple gene coding for the RdRp gene.
      if (not insert_result) {

        ExecEnv::log().warn("GeneFeature::getCodingSequences; Gene: {}, Duplicate coding feature: {}",
                            gene_ptr->id(),
                            leaf->featureText());
        gene_ptr->recusivelyPrintsubfeatures();
        return true;

      }

    }

  }

  if (not parent_feature_map.empty()) {

    std::shared_ptr<const TranscriptionSequence> coding_sequence(std::make_shared<const TranscriptionSequence>(gene_ptr,
                                                                                                               parent_ptr,
                                                                                                               parent_feature_map));

    result = result and sequence_array.insertSequence(coding_sequence);

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
