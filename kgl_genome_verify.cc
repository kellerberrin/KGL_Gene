//
// Created by kellerberrin on 12/11/17.
//


#include "kgl_patterns.h"
#include "kgl_genome_feature.h"
#include "kgl_genome_db.h"

namespace kgl = kellerberrin::genome;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ContigFeatures members.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////



void kgl::ContigFeatures::verifyFeatureHierarchy() {

  verifyContigOverlap();
  verifySubFeatureSuperFeatureDimensions();
  removeSubFeatureDuplicates();
  removeSuperFeatureDuplicates();
  verifySubFeatureDuplicates();
  verifySuperFeatureDuplicates();
  verifyCDSPhasePeptide();
  createGeneMap();

}

void kgl::ContigFeatures::createGeneMap() {

  // Clear the lookup table.
  gene_map_.clear();

  // Iterate through all the features looking for Gene features.
  for(const auto& feature : offset_feature_map_) {

    if(feature.second->isGene()) {

      ContigOffset_t end_offset = feature.second->sequence().end();
      gene_map_.insert(std::make_pair(end_offset, std::static_pointer_cast<GeneFeature>(feature.second)));

    }

  }

}


void kgl::ContigFeatures::verifyContigOverlap() {

  // If feature dimensions are [1, size] instead of [0, size) then assume that conversion from the
  // Gff convention of [1, size] has not been performed correctly during feature read from disk.
  // Adjust to [0, size) here.
  // Note that this suggests a problem with the (3rd party) Gff read functionality and should be addressed there.

  ContigSize_t contig_size = contigSize();

  for (auto feature_pair : id_feature_map_) {

    Feature &feature = *feature_pair.second;
    // Error if feature overlaps the and of the contig.
    // If [1,contig_size] then adjust to [0, contig_size)

    if (feature.sequence().begin() == 1) { // adjust to [0, size)

      FeatureSequence adj_sequence = feature.sequence();
      adj_sequence.begin(0);
      feature.sequence(adj_sequence);
      ExecEnv::log().warn("Contig: {} 1-offset features [1, {}], adjusted to zero-offset [0, {})",
                          contigId(), contig_size, contig_size);

    } else if (feature.sequence().end() > contig_size) { // No features larger than the contig.

      FeatureSequence adj_sequence = feature.sequence();
      adj_sequence.end(contig_size);
      feature.sequence(adj_sequence);
      ExecEnv::log().warn("Feature: {} [{}, {}) exceeds contig size :{} adjusted to [{}, {})",
                          feature.id(), feature.sequence().begin(), feature.sequence().end(),
                          contig_size, feature.sequence().begin(), contig_size);

    }

  } // for contig

}


void kgl::ContigFeatures::verifySubFeatureSuperFeatureDimensions() {


  // Check that sub-features fit within feature.
  for (auto feature_pair : id_feature_map_) {
    Feature &feature = *feature_pair.second;

    for (auto sub_feature_pair : feature.subFeatures()) {
      Feature &sub_feature = *sub_feature_pair.second;
      if (sub_feature.sequence().begin() < feature.sequence().begin()
          or sub_feature.sequence().end() > feature.sequence().end()) {

        kgl::ExecEnv::log().error("SubFeature: {}; [{}:{}] overlaps Feature {}; [{}:{}]",
                                  sub_feature.id(),
                                  sub_feature.sequence().begin(),
                                  sub_feature.sequence().end(),
                                  feature.id(),
                                  feature.sequence().begin(),
                                  feature.sequence().end());


      } // if sub-feature overlaps feature

    } // for all sub-features.

    // Check that features fit within super-feature.
    for (auto super_feature_pair : feature.superFeatures()) {
      Feature &super_feature = *super_feature_pair.second;
      if (feature.sequence().begin() < super_feature.sequence().begin()
          or feature.sequence().end() > super_feature.sequence().end()) {

        kgl::ExecEnv::log().error("Feature: {}; [{}:{}] overlaps SuperFeature {}; [{}:{}]",
                                  feature.id(),
                                  feature.sequence().begin(),
                                  feature.sequence().end(),
                                  super_feature.id(),
                                  super_feature.sequence().begin(),
                                  super_feature.sequence().end());

      } // If feature overlaps super-feature.

    } // for all super features.

  } // for all features.

}



void kgl::ContigFeatures::verifySubFeatureDuplicates() {

  for (auto feature_pair : id_feature_map_) {
    Feature &feature = *feature_pair.second;

    long duplicates = checkDuplicates(feature.subFeatures());

    if (duplicates > 0) {

      kgl::ExecEnv::log().warn("Feature: {}; has {} duplicate sub-features", feature.id(), duplicates);

    }

  }

}


void kgl::ContigFeatures::verifySuperFeatureDuplicates() {

  for (auto feature_pair : id_feature_map_) {
    Feature &feature = *feature_pair.second;

    long duplicates = checkDuplicates(feature.superFeatures());

    if (duplicates > 0) {

      kgl::ExecEnv::log().warn("Feature: {}; has {} duplicate super-features", feature.id(), duplicates);

    }

  }

}


void kgl::ContigFeatures::verifyCDSPhasePeptide() {

  // Iterate through all the features looking for Genes.
  size_t gene_count = 0;
  size_t ill_formed_genes = 0;
  size_t empty_genes = 0;

  ExecEnv::log().info("Verifying {} Gene structure using: {}",
                      contigId(), coding_sequence_.translationTableName());

  for(const auto& feature : offset_feature_map_) {

    if(feature.second->isGene()) {

      ++gene_count;
      const std::shared_ptr<const GeneFeature> gene_ptr = std::static_pointer_cast<const GeneFeature>(feature.second);
      const std::shared_ptr<const CodingSequenceArray> coding_seq_ptr = kgl::GeneFeature::getCodingSequences(gene_ptr);
      if (coding_seq_ptr->size() == 0) { // No CDS coding sequence available, try the EXON coding sequence

        ++empty_genes;

      }

      if (coding_seq_ptr->size() > 0) {

        if (not gene_ptr->verifyCDSPhase(coding_seq_ptr)) {

          ExecEnv::log().warn("Gene : {} Offset: {} Problem verifying CDS structure - gene sub-features print out",
                              gene_ptr->id(),
                              gene_ptr->sequence().begin());
          feature.second->recusivelyPrintsubfeatures();

        }

        if (not verifyCodingSequences(coding_seq_ptr)) {

          ++ill_formed_genes;

        }

      }

    }

  }

  ExecEnv::log().info("Verified {} found: {} Genes; Malformed Genes: {}", contigId(), gene_count, ill_formed_genes);

}


bool kgl::ContigFeatures::verifyCodingSequences(const std::shared_ptr<const CodingSequenceArray> coding_seq_ptr) const {

  bool result = true;

  if (coding_seq_ptr->size() == 0) {

    ExecEnv::log().error("codingSequence(), empty CodingSequenceArray");

  }

  for (const auto& sequence : coding_seq_ptr->getMap()) {

    if (sequence.second->getSortedCDS().empty()) {

      ExecEnv::log().error("codingSequence(), no corresponding CDS features found");
      continue;

    }

    std::shared_ptr<DNA5Sequence> coding_sequence_ptr = sequence_ptr_->codingSequence(sequence.second);

    if (not coding_sequence_.checkStartCodon(coding_sequence_ptr)) {

      ExecEnv::log().vwarn("No START codon Gene: {}, CDS parent (mRNA): {} | first codon: {}{}{}",
                           sequence.second->getGene()->id(),
                           sequence.second->getCDSParent()->id(),
                           coding_sequence_.firstCodon(coding_sequence_ptr).bases[0],
                           coding_sequence_.firstCodon(coding_sequence_ptr).bases[1],
                           coding_sequence_.firstCodon(coding_sequence_ptr).bases[2]);
//      gene_ptr->recusivelyPrintsubfeatures();
//      gene_ptr->printCDSvector(sorted_cds_vec);
      result = false;
    }
    if (not coding_sequence_.checkStopCodon(coding_sequence_ptr)) {

      ExecEnv::log().vwarn("No STOP codon: {} Gene: {}, CDS parent (mRNA): {} | last codon: {}{}{}",
                           (coding_sequence_.codonLength(coding_sequence_ptr)-1),
                           sequence.second->getGene()->id(),
                           sequence.second->getCDSParent()->id(),
                           coding_sequence_.lastCodon(coding_sequence_ptr).bases[0],
                           coding_sequence_.lastCodon(coding_sequence_ptr).bases[1],
                           coding_sequence_.lastCodon(coding_sequence_ptr).bases[2]);
//      gene_ptr->recusivelyPrintsubfeatures();
//      gene_ptr->printCDSvector(sorted_cds_vec);
      result = false;
    }
    size_t nonsense_index = coding_sequence_.checkNonsenseMutation(coding_sequence_ptr);
    if (nonsense_index > 0) {

      ExecEnv::log().vwarn("NONSENSE mutation codon:{} Gene: {}, CDS Parent (mRNA): {} | stop codon: {}{}{}",
                           nonsense_index,
                           sequence.second->getGene()->id(),
                           sequence.second->getCDSParent()->id(),
                           coding_sequence_.getCodon(coding_sequence_ptr, nonsense_index).bases[0],
                           coding_sequence_.getCodon(coding_sequence_ptr, nonsense_index).bases[1],
                           coding_sequence_.getCodon(coding_sequence_ptr, nonsense_index).bases[2]);
//      gene_ptr->recusivelyPrintsubfeatures();
//      gene_ptr->printCDSvector(sorted_cds_vec);
      result = false;
    }

  } // for cds group

  return result;

}


void kgl::ContigFeatures::removeSubFeatureDuplicates() {

  long duplicates_removed = 0;

  for (auto feature_pair : id_feature_map_) {
    Feature &feature = *feature_pair.second;

    duplicates_removed += deleteIterableDuplicates(feature.subFeatures());

  }

  if (duplicates_removed > 0) {

    kgl::ExecEnv::log().info("{} duplicate sub-features removed from contig: {}", duplicates_removed, contigId());

  }

}

void kgl::ContigFeatures::removeSuperFeatureDuplicates() {

  long duplicates_removed = 0;

  for (auto feature_pair : id_feature_map_) {
    Feature &feature = *feature_pair.second;

    duplicates_removed += deleteIterableDuplicates(feature.superFeatures());

  }

  if (duplicates_removed > 0) {

    kgl::ExecEnv::log().info("{} duplicate super-features removed from contig: {}", duplicates_removed, contigId());

  }

}

