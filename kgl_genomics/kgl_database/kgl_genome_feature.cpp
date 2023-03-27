//
// Created by kellerberrin on 10/10/17.
//

#include "kel_exec_env.h"
#include "kel_patterns.h"
#include "kgl_genome_feature.h"
#include "kgl_genome_contig.h"

#include <sstream>

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
     << contig()->contigId() << delimiter
     << "Feature Id:" << delimiter
     << id() << delimiter
     << "Type:" << delimiter
     << type() << delimiter
     << "Length:" << delimiter
     << sequence().length() << delimiter
     << "Offset:[" << delimiter
     << sequence().begin() << delimiter
     << "," << delimiter
     << sequence().end() << delimiter
     << ") Strand:" << delimiter
     << sequence().strandText() << delimiter
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

std::shared_ptr<const kgl::CodingSequenceArray>
kgl::GeneFeature::getCodingSequences(const std::shared_ptr<const GeneFeature>& gene) {

  if (gene->proteinCoding()) {

    std::shared_ptr<CodingSequenceArray> sequence_array_ptr(std::make_shared<CodingSequenceArray>(CodingSequenceType::PROTEIN));
    if (not getCodingSequences(Feature::CDS_TYPE_, gene, gene, sequence_array_ptr)) {

      ExecEnv::log().error("GeneFeature::getCodingSequences; Unable to retrieve CDS coding sequences for Gene: {}, type: {}", gene->id(), gene->type());

    }
    return sequence_array_ptr;

  } else if (gene->ncRNACoding()) {

    std::shared_ptr<CodingSequenceArray> sequence_array_ptr(std::make_shared<CodingSequenceArray>(CodingSequenceType::NCRNA));
    if (not getCodingSequences(Feature::MRNA_TYPE_, gene, gene, sequence_array_ptr)) {

      ExecEnv::log().error("GeneFeature::getCodingSequences; Unable to retrieve CDS coding sequences for Gene: {}, type: {}", gene->id(), gene->type());

    }
    return sequence_array_ptr;

  } else {

    std::shared_ptr<CodingSequenceArray> sequence_array_ptr(std::make_shared<CodingSequenceArray>(CodingSequenceType::PROTEIN));
    ExecEnv::log().error("GeneFeature::getCodingSequences; Gene: {}, Unknown gene type: {}", gene->id(), gene->type());
    return sequence_array_ptr;

  }

}

// This routine is recursive. Assumes all the CDS/EXONS are on the same sub-feature level.
bool kgl::GeneFeature::getCodingSequences(const std::string& coding_feature_type,
                                          const std::shared_ptr<const GeneFeature>& gene_ptr,
                                          const std::shared_ptr<const Feature>& coding_feature_parent_ptr,
                                          std::shared_ptr<CodingSequenceArray>& sequence_array_ptr) {

  bool result = true;
  TranscribedFeatureMap parent_cds;

  // loop through all sub_features
  for (const auto& [feature_id, feature_ptr] : coding_feature_parent_ptr->subFeatures()) {

    if (feature_ptr->type() == coding_feature_type) {

      auto [iter, insert_result] = parent_cds.insert(std::make_pair(feature_ptr->sequence().begin(), feature_ptr));

      // Some GFF files may have multiple coding features at the same logical level and the same begin offset.
      // This is true of the GFF supplied by NCBI for the SARS-COV-2 organism with multiple gene coding for the RdRp gene.
      if (not insert_result) {

        ExecEnv::log().warn("Gene: {}, Duplicate coding feature: {} at contig offset: {}",
                            gene_ptr->id(),
                            feature_id,
                            feature_ptr->sequence().begin());
//          gene_ptr->recusivelyPrintsubfeatures();


        // Call the coding sequence function recursively with the feature_ident argument set.
        result = result and getCodingSequences(coding_feature_type, gene_ptr, coding_feature_parent_ptr, sequence_array_ptr);

      }

    } else { // Assume feature is a higher feature such as mRNA and recursively call this function for sub-features.

      result = result and getCodingSequences(coding_feature_type, gene_ptr, feature_ptr, sequence_array_ptr);

    }

  } // for loop

  if (not parent_cds.empty()) {

    std::shared_ptr<const CodingSequence> coding_sequence(std::make_shared<const CodingSequence>(gene_ptr,
                                                                                                 coding_feature_parent_ptr,
                                                                                                 parent_cds));

    result = result and sequence_array_ptr->insertCodingSequence(coding_sequence);

  }

  return result;

}

