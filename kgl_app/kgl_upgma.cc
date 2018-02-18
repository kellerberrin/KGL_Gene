//
// Created by kellerberrin on 30/01/18.
//


#include <iomanip>
#include "kgl_upgma.h"


namespace kgl = kellerberrin::genome;




////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Mutates the contigs and then compares the mutated contigs using the Myer Hirschberg sequence comparison
// algorthim. This is linear in space - but quadratic in time. Need to find a faster comparison algorithm.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


kgl::DistanceType_t kgl::UPGMAContigDistance::distance(std::shared_ptr<const UPGMADistanceNode>  distance_node) const {

  std::shared_ptr<const UPGMAContigDistance> node_ptr = std::dynamic_pointer_cast<const UPGMAContigDistance>(distance_node);

  if (not node_ptr) {

    ExecEnv::log().error("distance(), Unexpected error, could not up-cast node pointer");
    return 1.0;

  }

  DistanceType_t total_distance = 0;
  for (auto contig : mutated_contigs_) {

    auto result = node_ptr->mutated_contigs_.find(contig.first);

    if (result != node_ptr->mutated_contigs_.end()) {

      ExecEnv::log().info("distance(), Comparing Genome: {}, Contig: {} with Genome: {}",
                          genome_variant_ptr_->genomeId(), contig.first, node_ptr->genome_variant_ptr_->genomeId());
      CompareDistance_t contig_distance = sequence_distance_->distance(contig.second, result->second);
      total_distance += static_cast<DistanceType_t>(contig_distance);
      ExecEnv::log().info("distance(), {} distance: {}", sequence_distance_->distanceType(), contig_distance);

    } else {

      ExecEnv::log().error("distance(), Unexpected error, could not find contig: {}", contig.first);

    }

  }

  return total_distance;

}



void kgl::UPGMAContigDistance::mutateContigs() {

  std::shared_ptr<const DNA5SequenceContig> reference_contig_ptr;
  std::shared_ptr<DNA5SequenceContig> mutant_contig_ptr;
  for (auto contig : genome_db_ptr_->getMap()) {

    if (not genome_variant_ptr_->mutantContig(contig.first, genome_db_ptr_, reference_contig_ptr, mutant_contig_ptr)) {

      ExecEnv::log().error("Unexpected error mutating contig; genome: {} , contig: {}", genome_variant_ptr_->genomeId(), contig.first);
      // Fail gracefully, just insert the reference contig.
      mutated_contigs_.insert(std::pair<ContigId_t, std::shared_ptr<const DNA5SequenceContig>>(contig.first, contig.second->sequence_ptr()));

    } else {

      mutated_contigs_.insert(std::pair<ContigId_t, std::shared_ptr<const DNA5SequenceContig>>(contig.first, mutant_contig_ptr));

    }

  }

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Mutates genome proteins and then compares the mutated proteins.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////



kgl::DistanceType_t kgl::UPGMAProteinDistance::distance(std::shared_ptr<const UPGMADistanceNode>  distance_node) const {

  std::shared_ptr<const UPGMAProteinDistance> node_ptr = std::dynamic_pointer_cast<const UPGMAProteinDistance>(distance_node);

  if (not node_ptr) {

    ExecEnv::log().error("distance(), Unexpected error, could not down-cast node pointer to UPGMAProteinDistance");
    return 1.0;

  }

  CompareDistance_t total_distance = 0;
  size_t gene_count = 0;

  for (auto protein : getMap()) {

    auto result = node_ptr->getMap().find(protein.first);

    if (result != node_ptr->getMap().end()) {

      CompareDistance_t contig_distance = sequence_distance_->distance(protein.second, result->second);
      total_distance += contig_distance;
      ++gene_count;

    } else {

      ExecEnv::log().error("distance(), Unexpected error, could not find Genome: {}, sequence: {}",
                           node_ptr->genomeId(), protein.first);

    }

  }

  ExecEnv::log().info("distance(), Genome: {}, Genome: {}; {} distance: {}, Gene Family: {}, Gene Count: {}",
                      genomeId(), node_ptr->genomeId(), sequence_distance_->distanceType(), total_distance, protein_family_, gene_count);
  return total_distance;

}


void kgl::UPGMAProteinDistance::mutateProteins() {

  mutated_proteins_.clear();
  for (auto contig : genome_db_ptr_->getMap()) {

    for (auto gene : contig.second->getGeneMap()) {

      if (protein_family_ != PROTEIN_FAMILY_WILDCARD) {

        std::shared_ptr<const OntologyRecord> ontology_record_ptr;
        if (ontology().getGafFeatureVector(gene.second->id(), ontology_record_ptr)) {

          if (ontology_record_ptr->symbolicReference() == protein_family_) {

            getProtein(gene.second);

          }
        }
      }
      else {

        getProtein(gene.second);

      }

    }

  }

}


void kgl::UPGMAProteinDistance::getProtein(std::shared_ptr<const GeneFeature> gene_ptr) {

  const std::shared_ptr<const CodingSequenceArray> coding_seq_ptr = kgl::GeneFeature::getCodingSequences(gene_ptr);
  for (auto sequence : coding_seq_ptr->getMap()) {

    std::shared_ptr<const ContigFeatures> contig_ptr = sequence.second->getGene()->contig();
    std::string gene_id = sequence.second->getGene()->id();
    std::string sequence_id = sequence.second->getCDSParent()->id();
    std::vector<std::shared_ptr<AminoSequence>> mutant_sequence_vector;
    std::shared_ptr<AminoSequence> reference_sequence;
    OffsetVariantMap variant_map;

    if (genome_variant_ptr_->mutantProteins(contig_ptr->contigId(),
                                             gene_id,
                                             sequence_id,
                                             genome_db_ptr_,
                                             variant_map,
                                             reference_sequence,
                                             mutant_sequence_vector)) {

      if (mutant_sequence_vector.empty()) {

        ExecEnv::log().warn("mutateProteins(), No mutant proteins generated for : genome: {} sequence: {}",
                            genome_variant_ptr_->genomeId(), sequence_id);

      } else if (mutant_sequence_vector.size() == 1) {

        mutated_proteins_[gene_id] = mutant_sequence_vector.front();

      } else {

        ExecEnv::log().warn("mutateProteins(), {} mutant proteins generated for : genome: {} sequence: {} only the first mutant is analyzed",
                            mutant_sequence_vector.size(), genome_variant_ptr_->genomeId(), sequence_id);

        mutated_proteins_[gene_id] = mutant_sequence_vector.front();

      }

    }

  }

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Mutates a single gene and compares to other selected genes
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool kgl::UPGMAGeneDistance::geneFamily(std::shared_ptr<const GeneFeature> gene_ptr,
                                        std::shared_ptr<const GenomeDatabase> genome_db_ptr,
                                        const std::string& protein_family) {

  std::shared_ptr<const OntologyRecord> ontology_record_ptr;
  if (genome_db_ptr->geneOntology().getGafFeatureVector(gene_ptr->id(), ontology_record_ptr)) {

    if (ontology_record_ptr->symbolicReference() == protein_family) {

      return true;

    }

  }

  return false;

}


void kgl::UPGMAGeneDistance::mutateProtein() {

  const std::shared_ptr<const CodingSequenceArray> coding_seq_ptr = kgl::GeneFeature::getCodingSequences(gene_ptr_);

  if (coding_seq_ptr->size() == 0) {

    ExecEnv::log().critical("mutateProtein(), Gene contains no coding sequence : genome: {} gene: {}",
                            genome_variant_ptr_->genomeId(), gene_ptr_->id());

  }

  std::shared_ptr<const CodingSequence> sequence = coding_seq_ptr->getFirst();
  std::shared_ptr<const ContigFeatures> contig_ptr = sequence->getGene()->contig();
  std::string gene_id = sequence->getGene()->id();
  std::string sequence_id = sequence->getCDSParent()->id();

  if (coding_seq_ptr->size() > 1) {

    ExecEnv::log().warn("mutateProtein(),  Genome: {} gene: {} contains: {} sequences using sequence: {}",
                        genome_variant_ptr_->genomeId(), gene_ptr_->id(), coding_seq_ptr->size(), sequence_id);

  }


  std::vector<std::shared_ptr<AminoSequence>> mutant_sequence_vector;
  std::shared_ptr<AminoSequence> reference_sequence;
  OffsetVariantMap variant_map;

  if (genome_variant_ptr_->mutantProteins(contig_ptr->contigId(),
                                           gene_id,
                                           sequence_id,
                                           genome_db_ptr_,
                                           variant_map,
                                           reference_sequence,
                                           mutant_sequence_vector)) {

    if (mutant_sequence_vector.empty()) {

      ExecEnv::log().error("mutateProtein(), No mutant proteins generated for : genome: {} sequence: {}",
                          genome_variant_ptr_->genomeId(), sequence_id);
      mutated_protein_ = reference_sequence;

    } else if (mutant_sequence_vector.size() == 1) {

      mutated_protein_ = mutant_sequence_vector.front();

    } else {

      ExecEnv::log().warn("mutateProtein(), {} mutant proteins generated for : genome: {} sequence: {} only the first mutant is analyzed",
                          mutant_sequence_vector.size(), genome_variant_ptr_->genomeId(), sequence_id);

      mutated_protein_ = mutant_sequence_vector.front();

    }

  }

}

kgl::DistanceType_t kgl::UPGMAGeneDistance::distance(std::shared_ptr<const UPGMADistanceNode>  distance_node) const {

  std::shared_ptr<const UPGMAGeneDistance> node_ptr = std::dynamic_pointer_cast<const UPGMAGeneDistance>(distance_node);

  if (not node_ptr) {

    ExecEnv::log().error("distance(), Unexpected error, could not down-cast node pointer to UPGMAGeneDistance");
    return 1.0;

  }


  CompareDistance_t contig_score = sequence_distance_->distance(mutated_protein_, node_ptr->mutated_protein_);
  DistanceType_t total_distance = static_cast<DistanceType_t>(contig_score);

  ExecEnv::log().info("distance(), Genome: {}, Gene: {}, Gene: {}; {} distance: {}, Gene Family: {}",
                      genome_variant_ptr_->genomeId(), gene_ptr_->id(), node_ptr->gene_ptr_->id(),
                      sequence_distance_->distanceType(), total_distance, protein_family_);

  return total_distance;

}


void kgl::UPGMAGeneDistance::writeNode(std::ofstream& outfile) const {

  std::stringstream ss;

  double contig_proportion = static_cast<double>(gene_ptr_->sequence().begin()) / static_cast<double>(gene_ptr_->contig()->contigSize());
  contig_proportion = contig_proportion * 100.0;
  ss << gene_ptr_->id() << "_" << std::setprecision(3) << contig_proportion;

  std::shared_ptr<const OntologyRecord> ontology_record_ptr;
  if (genome_db_ptr_->geneOntology().getGafFeatureVector(gene_ptr_->id(), ontology_record_ptr)) {

    ss << "_" << ontology_record_ptr->altSymbolicReference();

  }

  outfile << ss.str();

}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Compares a single gene between isolates.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void kgl::UPGMAATP4Distance::writeNode(std::ofstream& outfile) const {

  const std::shared_ptr<const CodingSequenceArray> coding_seq_ptr = kgl::GeneFeature::getCodingSequences(gene_ptr_);

  if (coding_seq_ptr->size() == 0) {

    ExecEnv::log().critical("write_node(), Gene contains no coding sequence : genome: {} gene: {}",
                            genome_variant_ptr_->genomeId(), gene_ptr_->id());

  }

  std::shared_ptr<const CodingSequence> sequence = coding_seq_ptr->getFirst();
  std::shared_ptr<const ContigFeatures> contig_ptr = sequence->getGene()->contig();
  std::string gene_id = sequence->getGene()->id();
  std::string sequence_id = sequence->getCDSParent()->id();

  if (coding_seq_ptr->size() > 1) {

    ExecEnv::log().warn("write_node(),  Genome: {} gene: {} contains: {} sequences using sequence: {}",
                        genome_variant_ptr_->genomeId(), gene_ptr_->id(), coding_seq_ptr->size(), sequence_id);

  }


  std::vector<std::shared_ptr<AminoSequence>> mutant_sequence_vector;
  std::shared_ptr<AminoSequence> reference_sequence;
  OffsetVariantMap variant_map;

  if (not genome_variant_ptr_->mutantProteins(contig_ptr->contigId(),
                                              gene_id,
                                              sequence_id,
                                              genome_db_ptr_,
                                              variant_map,
                                              reference_sequence,
                                              mutant_sequence_vector)) {

    ExecEnv::log().critical("write_node(), Cannot mutate sequence for : genome: {} gene: {}",
                            genome_variant_ptr_->genomeId(), gene_ptr_->id());

  }

  std::stringstream ss;

  ss << genome_variant_ptr_->genomeId();
  ss << "_" << SequenceManipulation().editItems(reference_sequence->getSequenceAsString(),
                                                mutated_protein_->getSequenceAsString(),
                                                '_',
                                                VariantOutputIndex::START_1_BASED);

  outfile << ss.str();

}


kgl::DistanceType_t kgl::UPGMAATP4Distance::zeroDistance(std::shared_ptr<const UPGMADistanceNode>  distance_node) const {

  std::shared_ptr<const UPGMAATP4Distance> node_ptr = std::dynamic_pointer_cast<const UPGMAATP4Distance>(distance_node);

  if (not node_ptr) {

    ExecEnv::log().error("distance(), Unexpected error, could not down-cast node pointer to UPGMAATP4Distance");
    return 1.0;

  }



  CompareDistance_t contig_score = zero_metric_ptr_->distance(mutated_protein_, node_ptr->mutated_protein_);
  auto total_distance = static_cast<DistanceType_t>(contig_score);

  ExecEnv::log().info("zero edit distance(), Genome: {}, Gene: {}, Gene: {}; {} distance: {}, Gene Family: {}",
                      genome_variant_ptr_->genomeId(), gene_ptr_->id(), node_ptr->gene_ptr_->id(),
                      zero_metric_ptr_->distanceType(), total_distance, protein_family_);

  return total_distance;

}
