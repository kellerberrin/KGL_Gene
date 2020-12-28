//
// Created by kellerberrin on 10/10/17.
//

#include "kel_exec_env.h"
#include "kel_patterns.h"
#include "kgl_genome_feature.h"
#include "kgl_genome_collection.h"

#include <sstream>

namespace kgl = kellerberrin::genome;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Feature members.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void kgl::Feature::addSubFeature(const FeatureIdent_t& sub_feature_id, std::shared_ptr<const Feature> sub_feature_ptr) {

  sub_features_.insert(std::make_pair(sub_feature_id, sub_feature_ptr));

}


void kgl::Feature::recusivelyPrintsubfeatures(long feature_level) const {

  for (const auto& feature : sub_features_) {

    ExecEnv::log().info("Level: {},  {}", feature_level, feature.second->featureText());

    feature.second->recusivelyPrintsubfeatures(feature_level + 1); // Recursive call for the sub-feature.

  }

}


std::string kgl::Feature::featureText() const {

  std::stringstream ss;

  ss << "Contig Id:"
     << contig()->contigId()
     << " Feature Id:"
     << id() 
     << " Type:"
     << featureType() 
     << " Length:"
     << sequence().length()
     << " Offset:["
     << sequence().begin() 
     << ", "
     << sequence().end() 
     << ") Strand:"
     << sequence().strandText();

  return ss.str();

}




std::shared_ptr<const kgl::Feature> kgl::Feature::getGene() const {

  // Recursively search upward

  if (hasSuperfeature()) {

    if (getSuperFeature()->isGene()) {

      return getSuperFeature();

    } else {

      return getSuperFeature()->getGene();

    }
  }

  return nullptr;

}


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Gene Feature members.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

std::shared_ptr<const kgl::CodingSequenceArray>
kgl::GeneFeature::getCodingSequences(std::shared_ptr<const GeneFeature> gene) {

  std::shared_ptr<CodingSequenceArray> sequence_array_ptr(std::make_shared<CodingSequenceArray>());
  if (not getCodingSequences(gene, gene, sequence_array_ptr)) {

    ExecEnv::log().error("getCodingSequences(), Unable to retrieve coding sequences for Gene: {}", gene->id());

  }

  return sequence_array_ptr;

}

// This routine is recursive. Assumes all the CDS/EXONS are on the same sub-feature level.
bool kgl::GeneFeature::getCodingSequences(std::shared_ptr<const GeneFeature> gene_ptr,
                                          std::shared_ptr<const Feature> cds_parent_ptr,
                                          std::shared_ptr<CodingSequenceArray>& sequence_array_ptr,
                                          const FeatureIdent_t& feature_ident /* optional argument*/ ) {

  bool result = true;
  SortedCDS parent_cds;

  // loop through all sub_features
  for (const auto& sub_feature : cds_parent_ptr->subFeatures()) {

    // if feature_ident is not set.
    if (feature_ident.empty()) {
 
      if (sub_feature.second->isCDS()) {


        auto insert = parent_cds.insert(std::make_pair(sub_feature.second->sequence().begin(),
                                        std::static_pointer_cast<const CDSFeature>(sub_feature.second)));

        // Some GFF files may have multiple coding features at the same logical level and the same begin offset.
        // The is true of the GFF supplied by NCBI for the SARS-COV-2 organism with multiple gene coding for the RdRp gene.

        if (not insert.second) {

          

          ExecEnv::log().warn("Gene: {}, Duplicate coding feature: {} at contig offset: {}",
                              gene_ptr->id(),
                              sub_feature.second->id(),
                              sub_feature.second->sequence().begin());
//          gene_ptr->recusivelyPrintsubfeatures();


          // Call the coding sequence function recursively with the feature_ident argument set.
          result = result and getCodingSequences(gene_ptr, cds_parent_ptr, sequence_array_ptr, sub_feature.second->id());

        }

      } else { // Assume feature is a higher feature such as mRNA and recursively call this function for sub-features.

        result = result and getCodingSequences(gene_ptr, sub_feature.second, sequence_array_ptr);

      }

    } else { // feature_ident specified.

      if (sub_feature.second->isCDS() and sub_feature.second->id() == feature_ident) {

        auto insert = parent_cds.insert(std::make_pair(sub_feature.second->sequence().begin(),
                                        std::static_pointer_cast<const CDSFeature>(sub_feature.second)));

        // Some GFF files may have multiple coding features at the same logical level and the same begin offset.
        // The is true of the GFF supplied by NCBI for the SARS-COV-2 organism with multiple gene coding for the RdRp gene.
        // todo: This logic should be changed so that the coding features are separated and become multiple gene sequences.
        // Just report a warning for now. 
        if (not insert.second) {

          ExecEnv::log().warn("Gene: {}, Duplicate coding feature: {} at contig offset: {}",
                              gene_ptr->id(),
                              sub_feature.second->id(),
                              sub_feature.second->sequence().begin());
//          gene_ptr->recusivelyPrintsubfeatures();
          result = false;

        }

      } 
      
    } // feature_ident specified.

  } // for loop

  if (not parent_cds.empty()) {

    std::shared_ptr<const CodingSequence> coding_sequence(std::make_shared<const CodingSequence>(gene_ptr,
                                                                                                 cds_parent_ptr,
                                                                                                 parent_cds));

    result = result and sequence_array_ptr->insertCodingSequence(coding_sequence);

  }

  return result;

}

