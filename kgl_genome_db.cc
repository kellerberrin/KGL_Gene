//
// Created by kellerberrin on 7/10/17.
//

#include "kgl_exec_env.h"
#include "kgl_patterns.h"
#include "kgl_genome_feature.h"
#include "kgl_genome_db.h"

namespace kgl = kellerberrin::genome;


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// ContigFeatures members.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


bool kgl::ContigFeatures::addFeature(std::shared_ptr<kgl::Feature>& feature_ptr) {

  id_feature_map_.insert(std::make_pair(feature_ptr->id(), feature_ptr));

  offset_feature_map_.insert(std::make_pair(feature_ptr->sequence().begin(), feature_ptr));

  return true;

}

bool kgl::ContigFeatures::findFeatureId(kgl::FeatureIdent_t& feature_id,
                                      std::vector<std::shared_ptr<kgl::Feature>>& feature_ptr_vec) {

  auto iter_pair = id_feature_map_.equal_range(feature_id);

  feature_ptr_vec.clear();
  for (auto iter = iter_pair.first; iter != iter_pair.second; ++iter) {

    feature_ptr_vec.emplace_back(iter->second);

  }

  return not feature_ptr_vec.empty();

}



void kgl::ContigFeatures::setupFeatureHierarchy() {

  // Remove all hierarchies for all features.
  for (auto feature_pair : id_feature_map_) {
    Feature& feature = *feature_pair.second;
    // Remove feature hierarchy.
    feature.clearHierachy();
  }

  // Establish or re-establish the hierarchies for all features.
  for (auto feature_pair : id_feature_map_) {
    Feature& feature = *feature_pair.second;
    // For each feature lookup a list of super_features
    std::vector<FeatureIdent_t> super_features;
    feature.getAttributes().getSuperFeatureIds(super_features);
    // Add parent pointers for the child and child pointers for the super_features.
    for (auto super_feature_id : super_features) {

      std::vector<std::shared_ptr<kgl::Feature>> super_feature_ptr_vec;
      if (not findFeatureId(super_feature_id, super_feature_ptr_vec)) {

        // Error; could not find super feature.
        kgl::ExecEnv::log().error("Feature: {}; Super Feature: {} does not exist", feature.id(), super_feature_id);

      }
      if (super_feature_ptr_vec.size() > 1) {

        // Warning, more than 1 super feature.
        kgl::ExecEnv::log().warn("Super Feature id: {} returned : {} Super Features",
                                 super_feature_id, super_feature_ptr_vec.size());

      }
      for (auto& super_feature_ptr : super_feature_ptr_vec) {

        feature.addSuperFeature(super_feature_id, super_feature_ptr);
        super_feature_ptr->addSubFeature(feature.id(), feature_pair.second);

      } // For all super_features with same id.

    } // For all parent ids.

  } // For all features.

}


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


bool kgl::ContigFeatures::findGenes(ContigOffset_t offset, GeneVector &gene_ptr_vec) const {

  gene_ptr_vec.clear();
  auto lb_result = gene_map_.lower_bound(offset);

  if (lb_result == gene_map_.end()) {

    return false;

  }

  auto result = gene_map_.equal_range(lb_result->first);

  for (auto it = result.first; it != result.second; ++it) {

    if (offset >= it->second->sequence().begin()) {

      gene_ptr_vec.emplace_back(it->second);

    }

  }

  return not gene_ptr_vec.empty();

}


// Convenience routine for tagging SNPs.
bool kgl::ContigFeatures::SNPMutation(std::shared_ptr<const CodingSequence> coding_seq_ptr,
                                      ContigOffset_t contig_offset,
                                      typename NucleotideColumn_DNA5::NucleotideType reference_base,
                                      typename NucleotideColumn_DNA5::NucleotideType mutant_base,
                                      ContigOffset_t& codon_offset,
                                      typename AminoAcidTypes::AminoType& reference_amino,
                                      typename AminoAcidTypes::AminoType& mutant_amino) const {


  return coding_sequence_.SNPMutation(coding_seq_ptr,
                                      sequence_ptr_,
                                      contig_offset,
                                      reference_base,
                                      mutant_base,
                                      codon_offset,
                                      reference_amino,
                                      mutant_amino);

}

// Convenience routine for Amino sequences.
std::shared_ptr<kgl::AminoSequence>
kgl::ContigFeatures::getAminoSequence(std::shared_ptr<const CodingSequence> coding_seq_ptr) const {

  return coding_sequence_.getAminoSequence(coding_seq_ptr, sequence_ptr_);

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// GenomeDatabase members.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool kgl::GenomeDatabase::addContigSequence(const kgl::ContigId_t& contig_id,
                                            std::shared_ptr<kgl::DNA5Sequence> sequence_ptr) {

  using ContigPtr = std::shared_ptr<kgl::ContigFeatures>;
  ContigPtr contig_ptr(std::make_shared<kgl::ContigFeatures>(contig_id, sequence_ptr));

  auto result = genome_sequence_map_.insert(std::make_pair(contig_id, std::move(contig_ptr)));

  return result.second;

}

bool kgl::GenomeDatabase::getContigSequence(const kgl::ContigId_t& contig_id,
                                             std::shared_ptr<ContigFeatures>& contig_ptr) const {

  auto result_iter = genome_sequence_map_.find(contig_id);

  if (result_iter != genome_sequence_map_.end()) {

    contig_ptr = result_iter->second;
    return true;

  }

  return false;

}

void kgl::GenomeDatabase::createVerifyGenomeDatabase() {

  setupFeatureHierarchy();
  verifyFeatureHierarchy();

}

void kgl::GenomeDatabase::setupFeatureHierarchy() {

  for (auto contig_pair : genome_sequence_map_) {

    contig_pair.second->setupFeatureHierarchy();

  }

}


void kgl::GenomeDatabase::verifyFeatureHierarchy() {

  for (auto contig_pair : genome_sequence_map_) {

    contig_pair.second->verifyFeatureHierarchy();

  }

}

void kgl::GenomeDatabase::setTranslationTable(size_t table) {

  for (auto contig_pair : genome_sequence_map_) {

    contig_pair.second->setTranslationTable(table);

  }

}


void kgl::GenomeDatabase::registerContigData(std::shared_ptr<kgl::ContigCountData>& contig_data_ptr) const {
// Create data blocks for each contig in the genome database
  for (const auto &contig_pair : genome_sequence_map_) {

    if (not contig_data_ptr->insertContig(contig_pair.first, contig_pair.second->sequence().length())) {

      kgl::ExecEnv::log().error("ContigCountData; attempted to add duplicate contig; {}", contig_pair.first);

    }

  }

}